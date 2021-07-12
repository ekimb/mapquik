use std::collections::VecDeque;

const seq_nt4_table: [u8; 256] =
   [0, 1, 2, 3,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  3, 3, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  3, 3, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4];
// copy from http://www.cse.yorku.ca/~oz/hash.html:

pub fn hash(mut key: u64) -> u64 {
        key = !key + (key << 21);
        key = key ^ key >> 24;
        key = (key + (key << 3)) + (key << 8);
        key = key ^ key >> 14;
        key = (key + (key << 2)) + (key << 4);
        key = key ^ key >> 28;

        key = key + (key << 31);
        return key;
}

pub fn update_window(mut q: &mut VecDeque<u64>, mut q_pos: &mut VecDeque<usize>, mut q_min_val: u64, mut q_min_pos: i32, new_strobe_hashval: u64, i: usize, mut new_minimizer: bool) -> (u64, i32, bool) {
    q.pop_front();
    let mut popped_index = q_pos.pop_front();
    q.push_back(new_strobe_hashval);
    q_pos.push_back(i);
    if q_min_pos == popped_index.unwrap() as i32 {
        q_min_val = u64::max_value();
        q_min_pos = i as i32;
        for j in q.len()..0 {
            if (q[j] < q_min_val) {
                q_min_val = q[j];
                q_min_pos = q_pos[j] as i32;
                new_minimizer = true;
            }
        }
    }
    else if new_strobe_hashval < q_min_val { // the new value added to queue is the new minimum
        q_min_val = new_strobe_hashval;
        q_min_pos = i as i32;
        new_minimizer = true;
    }
    (q_min_val, q_min_pos, new_minimizer)
}
    

pub fn extract_smers(k: usize, s: usize, seq: &[u8]) -> (Vec<u64>, Vec<usize>) {
    let smask : u64 = ((1 as u64) << 2*s) - 1;
    let kmask : u64 = ((1 as u64) << 2*k) - 1;
    let t = 3;
    let mut hash_count = 0;
    let mut seq_hashes = Vec::new();
    let mut pos_to_seq_coord = Vec::new();
    let mut qs = VecDeque::<u64>::new();
    let mut qs_pos = VecDeque::<usize>::new();
    let mut seq_len = seq.len();
    let mut qs_size = 0;
    let mut qs_min_val = u64::max_value();
    let mut qs_min_pos : i32 = -1;
    let mut l = 0;
    let mut xk : [u64; 2] = [0; 2];
    let mut xs : [u64; 2] = [0; 2];
    let mut kshift : u64 = (k as u64 - 1) * 2;
    let mut sshift : u64 = (s as u64 - 1) * 2;
    for i in 0..seq.len() {
        let mut c = seq_nt4_table[seq[i] as usize];
        if c < 4 {
            xk[0] = (xk[0] << 2 | c as u64) & kmask;                  // forward strand
            xk[1] = xk[1] >> 2 | ((3 - c) as u64) << kshift;  // reverse strand
            xs[0] = (xs[0] << 2 | c as u64) & smask;                  // forward strand
            xs[1] = xs[1] >> 2 | ((3 - c) as u64) << sshift;  // reverse strand
            l += 1;
            if l >= s {
                let mut ys : u64 = match xs[0] < xs[1]{
                    true => xs[0],
                    false => xs[1]
                };
                let hash_s = hash(ys);
                if qs_size < k - s {
                    qs.push_back(hash_s);
                    qs_pos.push_back(i - s+ 1);
                    qs_size += 1;
                }
                else if qs_size == k - s {
                    qs.push_back(hash_s);
                    qs_pos.push_back(i - s + 1);
                    qs_size += 1;
                    for j in 0..qs_size {
                        if qs[j] < qs_min_val {
                            qs_min_val = qs[j];
                            qs_min_pos = qs_pos[j] as i32;
                        }
                    }
                    if qs_min_pos == qs_pos[t-1] as i32 {
                        let mut yk : u64 = match xk[0] < xk[1]{
                            true => xk[0],
                            false => xk[1]
                        };
                        let mut hash_k = hash(yk);
                        seq_hashes.push(hash_k);
                        pos_to_seq_coord.push(i - k + 1);
                        hash_count += 1;
                    }
                }
                else {
                    let mut new_minimizer = false;
                    let tuple = update_window(&mut qs, &mut qs_pos, qs_min_val, qs_min_pos, hash_s, i - s + 1, new_minimizer);
                    qs_min_val = tuple.0; qs_min_pos = tuple.1; new_minimizer = tuple.2;
                    if qs_min_pos == qs_pos[t-1] as i32 { // occurs at t:th position in k-mer
                        let mut yk : u64 = match xk[0] < xk[1] {
                            true => xk[0],
                            false => xk[1]
                        };
                        let mut hash_k = hash(yk);
                        seq_hashes.push(hash_k);
                        pos_to_seq_coord.push(i - k + 1);
                        hash_count += 1;
                    }
                }
            }
        } else {
            qs_min_val = u64::max_value();
            qs_min_pos = -1;
            l = 0; xs[0] = 0; xs[1] = 0; xk[0] = 0; xk[1] = 0;
            qs_size = 0;
            qs.clear();
            qs_pos.clear();

        }
    }
    return (seq_hashes, pos_to_seq_coord)

} 