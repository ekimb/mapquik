// report.rs
// Builds the JSON coverage report (`--json`) directly from the mapping results.
//
// Schema stability: the keys of schema version 1 (`analysis_type`,
// `input_metrics`, `contig_coverage` and their keys) must stay unchanged;
// new fields may only be added.
//
// Everything here is pure data handling (no mapping logic) so it can be unit
// tested without running the mapper.

use crate::get_reader;
use seq_io::BaseRecord;
use serde::Serialize;
use serde_json::{json, Map, Value};
use std::collections::HashMap;
use std::fs::File;
use std::io::{self, BufWriter, Write};
use std::path::Path;

/// One PAF line produced by mapquik.
///
/// Columns 10 and 11 are kept as mapquik writes them (non-standard: `score`
/// and the reference length), only the coordinate columns are half-open.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PafRecord {
    pub query_name: String,
    pub query_len: u64,
    pub query_start: u64,
    pub query_end: u64,
    pub strand: char,
    pub target_name: String,
    pub target_len: u64,
    pub target_start: u64,
    pub target_end: u64,
    /// Column 10 (non-standard: number of matching seeds).
    pub score: u64,
    /// Column 11 (non-standard: reference length, kept for compatibility).
    pub block_len: u64,
    pub mapq: u64,
}

impl PafRecord {
    /// Formats the record as a 12-column PAF line (without the trailing newline).
    pub fn to_line(&self) -> String {
        format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            self.query_name,
            self.query_len,
            self.query_start,
            self.query_end,
            self.strand,
            self.target_name,
            self.target_len,
            self.target_start,
            self.target_end,
            self.score,
            self.block_len,
            self.mapq
        )
    }

    /// Parses a PAF line (at least the 12 mandatory columns are required).
    ///
    /// Returns `None` for malformed lines. mapquik always writes exactly 12
    /// columns, so in practice this only matters for externally produced PAFs.
    /// Kept as the inverse of [`PafRecord::to_line`] (used by the tests).
    #[allow(dead_code)]
    pub fn parse(line: &str) -> Option<PafRecord> {
        let mut fields = line.trim_end_matches(['\n', '\r']).split('\t');
        let mut next = || fields.next();
        let query_name = next()?.to_string();
        let query_len = next()?.parse().ok()?;
        let query_start = next()?.parse().ok()?;
        let query_end = next()?.parse().ok()?;
        let strand = match next()? {
            "+" => '+',
            "-" => '-',
            _ => return None,
        };
        let target_name = next()?.to_string();
        let target_len = next()?.parse().ok()?;
        let target_start = next()?.parse().ok()?;
        let target_end = next()?.parse().ok()?;
        let score = next()?.parse().ok()?;
        let block_len = next()?.parse().ok()?;
        let mapq = next()?.parse().ok()?;
        Some(PafRecord {
            query_name,
            query_len,
            query_start,
            query_end,
            strand,
            target_name,
            target_len,
            target_start,
            target_end,
            score,
            block_len,
            mapq,
        })
    }
}

/// Sorts and merges half-open `[start, end)` intervals.
///
/// Overlapping and adjacent (touching) intervals are merged, empty intervals
/// (`start == end`) are dropped. The result is sorted by start coordinate.
pub fn merge_intervals(intervals: &[(u64, u64)]) -> Vec<(u64, u64)> {
    let mut sorted: Vec<(u64, u64)> = intervals.iter().copied().filter(|(s, e)| e > s).collect();
    sorted.sort_unstable();
    let mut merged: Vec<(u64, u64)> = Vec::with_capacity(sorted.len());
    for (start, end) in sorted {
        match merged.last_mut() {
            Some(last) if start <= last.1 => {
                if end > last.1 {
                    last.1 = end;
                }
            }
            _ => merged.push((start, end)),
        }
    }
    merged
}

/// Number of bases covered by the union of the given half-open intervals.
pub fn covered_bases(intervals: &[(u64, u64)]) -> u64 {
    merge_intervals(intervals)
        .iter()
        .map(|(start, end)| end - start)
        .sum()
}

/// `value / total`, or `0.0` when `total` is zero.
pub fn ratio(value: u64, total: u64) -> f64 {
    if total == 0 {
        0.0
    } else {
        value as f64 / total as f64
    }
}

/// Rounds to `digits` decimal places (like Python's `round(x, digits)` for
/// non-tie values).
pub fn round_to(value: f64, digits: i32) -> f64 {
    let factor = 10f64.powi(digits);
    (value * factor).round() / factor
}

/// `value / total` as a percentage, rounded to two decimals; `0.0` if empty.
pub fn percent(value: u64, total: u64) -> f64 {
    if total == 0 {
        0.0
    } else {
        round_to((value as f64 / total as f64) * 100.0, 2)
    }
}

/// Contig name and length, in FASTA record order.
pub type ContigInfo = (String, u64);

/// Reads a FASTA/FASTQ reference a second time (sequentially) to obtain the
/// contig names and lengths in *file* order.
///
/// This is needed because reference indices are handed out from an atomic
/// counter in worker threads, so `ref_map` is not ordered, and the parallel
/// reader does not guarantee that its callbacks observe records in order.
/// Lengths are counted as non-whitespace bases, matching `normalize_seq()`.
pub fn contigs_in_file_order(path: &Path, fasta: bool) -> io::Result<Vec<ContigInfo>> {
    let (buf, _) = get_reader(path);
    let mut contigs = Vec::new();
    if fasta {
        let mut reader = seq_io::fasta::Reader::new(buf);
        while let Some(record) = reader.next() {
            let record = record.map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
            let name = record.id().unwrap_or("").to_string();
            let len = record
                .seq_lines()
                .map(|line| line.iter().filter(|b| !b.is_ascii_whitespace()).count() as u64)
                .sum();
            contigs.push((name, len));
        }
    } else {
        let mut reader = seq_io::fastq::Reader::new(buf);
        while let Some(record) = reader.next() {
            let record = record.map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
            let name = record.id().unwrap_or("").to_string();
            let len = record
                .seq()
                .iter()
                .filter(|b| !b.is_ascii_whitespace())
                .count() as u64;
            contigs.push((name, len));
        }
    }
    Ok(contigs)
}

/// Per-contig accumulated coverage.
#[derive(Debug, Clone)]
pub struct ContigCoverage {
    pub name: String,
    pub length: u64,
    pub total_aligned_bases: u64,
    pub mapped_reads: u64,
    pub mapped_reads_mapq60: u64,
    /// Reference intervals `[target_start, target_end)` of the PAF lines.
    pub intervals: Vec<(u64, u64)>,
}

/// Accumulates the report data while the PAF is written.
///
/// Not thread safe on purpose: it is used from the main thread (seq_io path) or
/// from the single PAF writer thread (parallelfastx path). Input read counters
/// for the parallelfastx path are collected with atomics and merged in with
/// [`CoverageAccumulator::set_input_totals`].
#[derive(Debug)]
pub struct CoverageAccumulator {
    contigs: Vec<ContigCoverage>,
    index: HashMap<String, usize>,
    pub input_reads_total: u64,
    pub input_bases_total: u64,
    pub reads_mapped_count: u64,
    pub reads_mapped_bases: u64,
}

impl CoverageAccumulator {
    /// `contigs` must be in report (= reference FASTA) order.
    pub fn new(contigs: Vec<ContigInfo>) -> Self {
        let mut index = HashMap::with_capacity(contigs.len());
        let contigs = contigs
            .into_iter()
            .enumerate()
            .map(|(i, (name, length))| {
                index.insert(name.clone(), i);
                ContigCoverage {
                    name,
                    length,
                    total_aligned_bases: 0,
                    mapped_reads: 0,
                    mapped_reads_mapq60: 0,
                    intervals: Vec::new(),
                }
            })
            .collect();
        CoverageAccumulator {
            contigs,
            index,
            input_reads_total: 0,
            input_bases_total: 0,
            reads_mapped_count: 0,
            reads_mapped_bases: 0,
        }
    }

    /// Counts one query record (mapped or not).
    pub fn add_input_read(&mut self, len: u64) {
        self.input_reads_total += 1;
        self.input_bases_total += len;
    }

    /// Overwrites the input totals (used after the parallelfastx atomics).
    pub fn set_input_totals(&mut self, reads: u64, bases: u64) {
        self.input_reads_total = reads;
        self.input_bases_total = bases;
    }

    /// Counts one PAF line. mapquik emits at most one line per read.
    pub fn add_paf_record(&mut self, record: &PafRecord) {
        let idx = match self.index.get(&record.target_name) {
            Some(&idx) => idx,
            None => {
                // Should not happen for records produced by this run; keep the
                // data instead of dropping it.
                let idx = self.contigs.len();
                self.index.insert(record.target_name.clone(), idx);
                self.contigs.push(ContigCoverage {
                    name: record.target_name.clone(),
                    length: record.target_len,
                    total_aligned_bases: 0,
                    mapped_reads: 0,
                    mapped_reads_mapq60: 0,
                    intervals: Vec::new(),
                });
                idx
            }
        };
        let contig = &mut self.contigs[idx];
        contig.total_aligned_bases += record.target_end.saturating_sub(record.target_start);
        contig.mapped_reads += 1;
        if record.mapq == 60 {
            contig.mapped_reads_mapq60 += 1;
        }
        contig
            .intervals
            .push((record.target_start.min(record.target_end), record.target_end.max(record.target_start)));
        self.reads_mapped_count += 1;
        self.reads_mapped_bases += record.query_len;
    }

    /// Builds the JSON report.
    pub fn build(&self, meta: &ReportMeta) -> Value {
        let total_assembly_length_bp: u64 = self.contigs.iter().map(|c| c.length).sum();
        let total_aligned_bases: u64 = self.contigs.iter().map(|c| c.total_aligned_bases).sum();
        let total_covered_bases: u64 = self.contigs.iter().map(|c| covered_bases(&c.intervals)).sum();

        let reads_unmapped_count = self.input_reads_total.saturating_sub(self.reads_mapped_count);

        let input_metrics = json!({
            "input_reads_total": self.input_reads_total,
            "reads_mapped_count": self.reads_mapped_count,
            "reads_unmapped_count": reads_unmapped_count,
            "reads_unmapped_percent": percent(reads_unmapped_count, self.input_reads_total),
            "total_assembly_length_bp": total_assembly_length_bp,
            "reads_mapped_percent": percent(self.reads_mapped_count, self.input_reads_total),
            "input_bases_total": self.input_bases_total,
            "reads_mapped_bases": self.reads_mapped_bases,
            "total_aligned_bases": total_aligned_bases,
            "average_depth": round_to(ratio(total_aligned_bases, total_assembly_length_bp), 2),
            "covered_bases": total_covered_bases,
            "breadth_fraction": round_to(ratio(total_covered_bases, total_assembly_length_bp), 4),
        });

        let mut contig_coverage = Map::new();
        for contig in &self.contigs {
            let contig_covered_bases = covered_bases(&contig.intervals);
            contig_coverage.insert(
                contig.name.clone(),
                json!({
                    "length_bp": contig.length,
                    "total_aligned_bases": contig.total_aligned_bases,
                    "average_depth": round_to(ratio(contig.total_aligned_bases, contig.length), 2),
                    "covered_bases": contig_covered_bases,
                    "breadth_fraction": round_to(
                        ratio(contig_covered_bases, contig.length),
                        4
                    ),
                    "mapped_reads": contig.mapped_reads,
                    "mapped_reads_mapq60": contig.mapped_reads_mapq60,
                }),
            );
        }

        json!({
            "analysis_type": "Mapquik_Coverage_Report",
            "schema_version": 2,
            "mapquik": meta.to_json(),
            "input_metrics": input_metrics,
            "contig_coverage": Value::Object(contig_coverage),
        })
    }
}

/// Parameters recorded in the report, using the command line option names.
#[derive(Debug, Clone)]
pub struct Parameters {
    pub preset: String,
    pub k: usize,
    pub l: usize,
    pub density: f64,
    pub chain: usize,
    pub seed: usize,
    pub gap_diff: usize,
    pub b: usize,
    pub q: usize,
    pub threads: usize,
    pub low_memory: bool,
    pub hpc: bool,
    pub simd: bool,
    pub parallelfastx: bool,
    pub debug: bool,
}

impl Parameters {
    fn to_json(&self) -> Value {
        json!({
            "preset": self.preset,
            "k": self.k,
            "l": self.l,
            "density": self.density,
            "chain": self.chain,
            "seed": self.seed,
            "gap_diff": self.gap_diff,
            "b": self.b,
            "q": self.q,
            "threads": self.threads,
            "low_memory": self.low_memory,
            "hpc": self.hpc,
            "simd": self.simd,
            "parallelfastx": self.parallelfastx,
            "debug": self.debug,
        })
    }
}

/// Report provenance: version, inputs and parameters.
#[derive(Debug, Clone)]
pub struct ReportMeta {
    pub version: String,
    pub reads: String,
    pub reference: String,
    pub parameters: Parameters,
}

impl ReportMeta {
    fn to_json(&self) -> Value {
        json!({
            "version": self.version,
            "reads": self.reads,
            "reference": self.reference,
            "parameters": self.parameters.to_json(),
        })
    }
}

/// Writes the report as pretty-printed JSON with a 4-space indent (same layout
/// as Python's `json.dump(..., indent=4)`).
pub fn write_report(
    path: &Path,
    accumulator: &CoverageAccumulator,
    meta: &ReportMeta,
) -> io::Result<()> {
    let value = accumulator.build(meta);
    let file = File::create(path)?;
    let mut writer = BufWriter::new(file);
    let formatter = serde_json::ser::PrettyFormatter::with_indent(b"    ");
    let mut serializer = serde_json::Serializer::with_formatter(&mut writer, formatter);
    value
        .serialize(&mut serializer)
        .map_err(|e| io::Error::new(io::ErrorKind::Other, e))?;
    writer.flush()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn meta() -> ReportMeta {
        ReportMeta {
            version: "0.0.0".to_string(),
            reads: "reads.fa".to_string(),
            reference: "ref.fa".to_string(),
            parameters: Parameters {
                preset: "hifi".to_string(),
                k: 8,
                l: 16,
                density: 0.01,
                chain: 4,
                seed: 11,
                gap_diff: 100,
                b: 1,
                q: 200,
                threads: 8,
                low_memory: false,
                hpc: true,
                simd: true,
                parallelfastx: false,
                debug: false,
            },
        }
    }

    fn contigs() -> Vec<ContigInfo> {
        vec![("chr2".to_string(), 1000), ("chr1".to_string(), 500)]
    }

    fn record(query: &str, q_len: u64, target: &str, target_start: u64, target_end: u64, mapq: u64) -> PafRecord {
        PafRecord {
            query_name: query.to_string(),
            query_len: q_len,
            query_start: 0,
            query_end: q_len,
            strand: '+',
            target_name: target.to_string(),
            target_len: 0,
            target_start,
            target_end,
            score: 20,
            block_len: 0,
            mapq,
        }
    }

    #[test]
    fn merge_overlapping_intervals() {
        assert_eq!(
            merge_intervals(&[(10, 20), (15, 30), (0, 5)]),
            vec![(0, 5), (10, 30)]
        );
        assert_eq!(covered_bases(&[(10, 20), (15, 30), (0, 5)]), 25);
    }

    #[test]
    fn merge_adjacent_intervals() {
        assert_eq!(merge_intervals(&[(0, 10), (10, 20)]), vec![(0, 20)]);
        assert_eq!(covered_bases(&[(0, 10), (10, 20)]), 20);
        // Touching but unsorted.
        assert_eq!(merge_intervals(&[(20, 30), (0, 10), (10, 20)]), vec![(0, 30)]);
    }

    #[test]
    fn merge_nested_intervals() {
        assert_eq!(merge_intervals(&[(0, 100), (10, 20), (30, 40)]), vec![(0, 100)]);
        assert_eq!(covered_bases(&[(0, 100), (10, 20)]), 100);
    }

    #[test]
    fn merge_disjoint_intervals() {
        assert_eq!(
            merge_intervals(&[(100, 110), (0, 10), (50, 60)]),
            vec![(0, 10), (50, 60), (100, 110)]
        );
        assert_eq!(covered_bases(&[(100, 110), (0, 10), (50, 60)]), 30);
    }

    #[test]
    fn merge_ignores_empty_intervals() {
        assert_eq!(merge_intervals(&[(5, 5), (7, 7)]), Vec::new());
        assert_eq!(covered_bases(&[]), 0);
    }

    #[test]
    fn paf_roundtrip() {
        let line = "read1\t100\t0\t100\t+\tchr1\t500\t10\t110\t42\t500\t60";
        let record = PafRecord::parse(line).unwrap();
        assert_eq!(record.query_name, "read1");
        assert_eq!(record.target_name, "chr1");
        assert_eq!(record.mapq, 60);
        assert_eq!(record.to_line(), line);
        assert!(PafRecord::parse("read1\t100").is_none());
        assert!(PafRecord::parse("read1\t100\t0\t100\t*\tchr1\t500\t10\t110\t42\t500\t60").is_none());
    }

    #[test]
    fn report_has_wrapper_keys_and_extra_fields() {
        let mut accumulator = CoverageAccumulator::new(contigs());
        // Two mapped reads on chr1: one overlapping (nested), one adjacent.
        accumulator.add_input_read(100);
        accumulator.add_input_read(50);
        accumulator.add_input_read(30); // unmapped
        accumulator.add_paf_record(&record("read1", 100, "chr1", 0, 100, 60));
        accumulator.add_paf_record(&record("read2", 50, "chr1", 50, 120, 0));

        let report = accumulator.build(&meta());

        assert_eq!(report["analysis_type"], "Mapquik_Coverage_Report");
        assert_eq!(report["schema_version"], 2);
        assert_eq!(report["mapquik"]["version"], "0.0.0");
        assert_eq!(report["mapquik"]["parameters"]["k"], 8);
        assert_eq!(report["mapquik"]["parameters"]["preset"], "hifi");

        let metrics = &report["input_metrics"];
        assert_eq!(metrics["input_reads_total"], 3);
        assert_eq!(metrics["reads_mapped_count"], 2);
        assert_eq!(metrics["reads_unmapped_count"], 1);
        assert_eq!(metrics["reads_unmapped_percent"], 33.33);
        assert_eq!(metrics["reads_mapped_percent"], 66.67);
        assert_eq!(metrics["total_assembly_length_bp"], 1500);
        assert_eq!(metrics["input_bases_total"], 180);
        assert_eq!(metrics["reads_mapped_bases"], 150);
        assert_eq!(metrics["total_aligned_bases"], 170);
        assert_eq!(metrics["average_depth"], 0.11); // 170 / 1500 = 0.1133
        assert_eq!(metrics["covered_bases"], 120);
        assert_eq!(metrics["breadth_fraction"], 0.08); // 120 / 1500

        let chr1 = &report["contig_coverage"]["chr1"];
        assert_eq!(chr1["length_bp"], 500);
        assert_eq!(chr1["total_aligned_bases"], 170);
        assert_eq!(chr1["average_depth"], 0.34);
        assert_eq!(chr1["covered_bases"], 120);
        assert_eq!(chr1["breadth_fraction"], 0.24);
        assert_eq!(chr1["mapped_reads"], 2);
        assert_eq!(chr1["mapped_reads_mapq60"], 1);

        let chr2 = &report["contig_coverage"]["chr2"];
        assert_eq!(chr2["length_bp"], 1000);
        assert_eq!(chr2["total_aligned_bases"], 0);
        assert_eq!(chr2["average_depth"], 0.0);
        assert_eq!(chr2["covered_bases"], 0);
        assert_eq!(chr2["breadth_fraction"], 0.0);
        assert_eq!(chr2["mapped_reads"], 0);
        assert_eq!(chr2["mapped_reads_mapq60"], 0);

        // Reference order is kept, even though records arrived for chr1 first.
        let order: Vec<&String> = report["contig_coverage"]
            .as_object()
            .unwrap()
            .keys()
            .collect();
        assert_eq!(order, vec!["chr2", "chr1"]);
    }

    #[test]
    fn report_without_reads() {
        let accumulator = CoverageAccumulator::new(contigs());
        let report = accumulator.build(&meta());

        let metrics = &report["input_metrics"];
        assert_eq!(metrics["input_reads_total"], 0);
        assert_eq!(metrics["reads_mapped_count"], 0);
        assert_eq!(metrics["reads_unmapped_count"], 0);
        // No division by zero.
        assert_eq!(metrics["reads_unmapped_percent"], 0.0);
        assert_eq!(metrics["reads_mapped_percent"], 0.0);
        assert_eq!(metrics["average_depth"], 0.0);
        assert_eq!(metrics["breadth_fraction"], 0.0);
        assert_eq!(metrics["total_assembly_length_bp"], 1500);

        assert_eq!(report["contig_coverage"]["chr1"]["average_depth"], 0.0);
        assert_eq!(report["contig_coverage"]["chr1"]["breadth_fraction"], 0.0);
    }

    #[test]
    fn zero_length_contig_does_not_divide_by_zero() {
        let mut accumulator = CoverageAccumulator::new(vec![("empty".to_string(), 0)]);
        accumulator.add_paf_record(&record("read1", 10, "empty", 0, 0, 60));
        let report = accumulator.build(&meta());
        assert_eq!(report["contig_coverage"]["empty"]["average_depth"], 0.0);
        assert_eq!(report["contig_coverage"]["empty"]["breadth_fraction"], 0.0);
        assert_eq!(report["input_metrics"]["average_depth"], 0.0);
        assert_eq!(report["input_metrics"]["breadth_fraction"], 0.0);
    }

    #[test]
    fn rounding_matches_python_style() {
        assert_eq!(round_to(0.113333, 2), 0.11);
        assert_eq!(round_to(0.116666, 2), 0.12);
        assert_eq!(round_to(1.0 / 3.0, 4), 0.3333);
        assert_eq!(percent(1, 3), 33.33);
        assert_eq!(percent(2, 3), 66.67);
        assert_eq!(percent(1, 0), 0.0);
    }

    #[test]
    fn contigs_are_read_in_file_order() {
        let path = std::env::temp_dir().join(format!("mapquik_report_test_{}.fa", std::process::id()));
        std::fs::write(
            &path,
            b">zz first record\nACGT\nACGT\n>aa\nTTTT\n>mm\nGGGG\nCCCC\n",
        )
        .unwrap();
        let contigs = contigs_in_file_order(&path, true).unwrap();
        std::fs::remove_file(&path).ok();
        assert_eq!(
            contigs,
            vec![
                ("zz".to_string(), 8),
                ("aa".to_string(), 4),
                ("mm".to_string(), 8)
            ]
        );
    }

    #[test]
    fn contigs_from_fastq_reference() {
        let path = std::env::temp_dir().join(format!("mapquik_report_test_{}.fq", std::process::id()));
        std::fs::write(&path, b"@zz\nACGTACGT\n+\nIIIIIIII\n@aa\nTTTT\n+\nIIII\n").unwrap();
        let contigs = contigs_in_file_order(&path, false).unwrap();
        std::fs::remove_file(&path).ok();
        assert_eq!(contigs, vec![("zz".to_string(), 8), ("aa".to_string(), 4)]);
    }
}
