use std::str::FromStr;
use std::num::ParseIntError;
use regex::Regex;

// Define a tuple struct to hold a CIGAR operation (u32, char)
#[derive(Debug, PartialEq)]
struct CigarUnit(u32, char);

impl FromStr for CigarUnit {
    type Err = ParseIntError;
    
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let (num_str, op_str) = s.split_at(s.len() - 1);
        let num = num_str.parse::<u32>()?;
        let op = op_str.chars().next().unwrap();
        Ok(CigarUnit(num, op))
    }
}

pub fn merge_cigar_strings(cigar_strings: Vec<String>) -> String {
    let mut cigar_parts: Vec<(usize, char)> = Vec::new();

    for cigar in cigar_strings {
        let re = Regex::new(r"([0-9]+)([A-Za-z=])").unwrap();

        for cap in re.captures_iter(&cigar) {
            let num = cap[1].parse().unwrap();
            let ch = cap[2].chars().next().unwrap();
            if let Some((last_num, last_ch)) = cigar_parts.last_mut() {
                if *last_ch == ch {
                    *last_num += num;
                } else {
                    cigar_parts.push((num, ch));
                }
            } else {
                cigar_parts.push((num, ch));
            }
        }
    }

    let mut output = String::new();
    for (num, ch) in cigar_parts {
        output.push_str(&format!("{}{}", num, ch));
    }

    output
}
