


import os,sys
import argparse

import pysam


def read_sam(sam_file):
    SAM_file = pysam.AlignmentFile(sam_file, "r", check_sq=False)
    read_positions = {} # acc -> [ref_id, ref_start, refstop]


    for read in SAM_file.fetch(until_eof=True):
        if read.flag == 0 or read.flag == 16:
            # print(read.query_name, len(read_positions))
            read_positions[read.query_name] = (read.reference_name, read.reference_start, read.reference_end)
        elif read.flag == 4:
            read_positions[read.query_name] = False

    return read_positions     


def read_paf(paf_file):
    read_positions = {} # acc -> [ref_id, ref_start, refstop]
    mapped_to_multiple_pos = 0
    for line in open(paf_file, 'r'):
        vals = line.split()
        read_acc, ref_name, reference_start, reference_end  = vals[0], vals[5], int(vals[7]), int(vals[8])

        if read_acc in read_positions:
            mapped_to_multiple_pos += 1
            continue
        else:
            read_positions[read_acc] = (ref_name, reference_start, reference_end)
    return read_positions, mapped_to_multiple_pos


def overlap(q_a, q_b, p_a, p_b):
    assert q_a <= q_b and p_a <= p_b
    # if (q_a == q_b) or (p_a == p_b):
    #     print("Cigar bug")
    return  (p_a <= q_a <= p_b) or (p_a <= q_b <= p_b) or (q_a <= p_a <= q_b) or (q_a <= p_b <= q_b)

def get_stats(truth, predicted):

    # nr_aligned = len(predicted)
    nr_total = len(truth)
    unaligned = 0 
    nr_aligned = 0
    over_mapped = 0
    correct = 0
    for read_acc in predicted:
        if not truth[read_acc]:
            over_mapped += 1
            continue
        if not predicted[read_acc]:
            unaligned += 1
            continue

        nr_aligned += 1

        pred_ref_id, pred_start, pred_stop = predicted[read_acc]
        true_ref_id, true_start, true_stop = truth[read_acc]
        # print(read_acc, pred_start, pred_stop, true_start, true_stop)
        if pred_ref_id == true_ref_id and overlap(pred_start, pred_stop, true_start, true_stop):
            correct += 1
            # print(read_acc)
        else:
            pass
            # print(read_acc, pred_ref_id, pred_start, pred_stop, true_ref_id, true_start, true_stop )

    return 100*(nr_aligned/nr_total), 100*correct/nr_aligned, over_mapped


def main(args):

    truth = read_sam(args.truth)

    if args.predicted_sam:
        predicted = read_sam(args.predicted_sam)
    elif args.predicted_paf:
        predicted, mapped_to_multiple_pos = read_paf(args.predicted_paf)
        # print("Number of reads mapped to several positions (using first pos):", mapped_to_multiple_pos)

    percent_aligned, percent_correct, over_mapped = get_stats(truth, predicted)
    out_str = "{0},{1},{2}".format(round(percent_aligned, 3), round(percent_correct, 3), over_mapped)
    print(out_str)
    # print("Percentage aligned: {0}".format(round(percent_aligned, 3)))
    # print("Accuracy: {0}".format(round(percent_correct, 3)))
    # print("Over-aligned (unmapped in grough truth file): {0}".format(over_mapped))



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Calc identity", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--truth', type=str, default=False, help='True coordinates (SAM)')
    parser.add_argument('--predicted_sam', type=str, default="", help='Perdicted coordinates (SAM/PAF)')
    parser.add_argument('--predicted_paf', type=str, default="", help='Perdicted coordinates (SAM/PAF)')
    parser.add_argument('--outfile', type=str, default=None, help='Path to file.')
    # parser.set_defaults(which='main')
    args = parser.parse_args()



    if len(sys.argv)==1:
        parser.print_help()
        sys.exit()

    main(args)