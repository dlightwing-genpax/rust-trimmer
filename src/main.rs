use std::io::prelude::*;
use std::fs::File;
use std::time::*;
use std::cmp::*;
use std::path::Path;
use std::fs;
use flate2::read::MultiGzDecoder;
use flate2::write::GzEncoder;
use flate2::Compression;
use itertools::Itertools;
use rayon::prelude::*;
use bitvec::prelude::*;
use clap::Parser;

type BV = BitVec<u8, Msb0>;
type BS = BitSlice<u8, Msb0>;

fn render(mut k: usize, klen: usize) -> String {
    let mut bases = vec![0 as u8; klen];
    for i in 0..klen {
        bases[klen-i-1] = (k % 4) as u8;
        k >>= 2;
    }
    bases.iter().map(|b| 
        match b {0=>'A', 1=>'T', 2=>'C', 3=>'G', _=>'N'}
    ).collect::<String>()
}

fn unpack(slice: &BS) -> usize {
    if slice.len() == 0 {
        return 0;
    }
    slice.load_be::<usize>()
}

type ReadData = (BV, BV, BV, String, String, String);


fn load_readset(fq1: &str, fq2: &str) -> Vec<[ReadData; 2]> {
    let reads = [fq1, fq2].into_par_iter().map(|fq| {
        let mut f = File::open(fq).unwrap();
        let mut s_gz = Vec::new();
        f.read_to_end(&mut s_gz).unwrap();

        let mut decoder = MultiGzDecoder::new(s_gz.as_slice());
        let mut s = String::new();
        if decoder.read_to_string(&mut s).err().is_some() {
            panic!("{fq} is invalid gzip.");
        }

        s.lines().chunks(4)
            .into_iter().map(|c| c.collect::<Vec<_>>()).collect::<Vec<_>>()
            .into_par_iter().map(|read| {
                let l = read[1].len();
                let mut rtn = (
                    BV::repeat(false, l),
                    BV::repeat(false, l),
                    BV::repeat(false, l),
                    read[0].to_string(),
                    read[1].to_string(),
                    read[3].to_string(),
                );
            
                for (idx, (b, q)) in read[1].chars().map(|c| 
                    match c {'A'=>0, 'T'=>1, 'C'=>2, 'G'=>3, _=>4}
                ).into_iter().zip(read[3].chars()).enumerate() {
                    if b == 4 {
                        continue;
                    } else {
                        rtn.0.set(idx, true);
                        rtn.1.set(idx, b >> 1 == 1);
                        rtn.2.set(idx, b % 2 == 1);
                    }
                }
                rtn
            }).collect::<Vec<_>>()
    }).collect::<Vec<_>>();

    // put corresponding r1 & r2 together
    reads[0].par_iter().zip(reads[1].par_iter()).map(|(r1, r2)| 
        [r1.clone(), r2.clone()]
    ).collect::<Vec<_>>()
}


fn get_adapter_candidates(read: &[ReadData; 2], inverses: &Vec<usize>, klen: usize, max_diff: usize) 
    -> Vec<(bool, usize, usize)> {
    // outputs the first klen bases of all possible adapters as (side, bit1, bit2)
    let mut adapters = Vec::new();

    if read[0].0.len() < 32 || read[1].0.len() < 32 {
        return vec![];
    }

    for is in 0 .. max(read[0].0.len(), read[1].0.len()) + 1 - klen {
        // is => insert size

        if is + klen <= min(read[0].0.len(), read[1].0.len()) {
            let force_klen = min(is, klen);
            
            let q0 = unpack(&read[0].0[is-force_klen .. is]) & inverses[unpack(&read[1].0[0 .. force_klen])];
            let q1 = unpack(&read[0].1[is-force_klen .. is]) ^ inverses[unpack(&read[1].1[0 .. force_klen])];
            let q2 = unpack(&read[0].2[is-force_klen .. is]) ^ !inverses[unpack(&read[1].2[0 .. force_klen])];
            let qscore = usize::count_ones(q0 & !q1 & !q2) as usize;   
            let q_ok = read[0].0[is .. is+klen].count_ones() == klen;

            let r0 = unpack(&read[1].0[is-force_klen .. is]) & inverses[unpack(&read[0].0[0 .. force_klen])];
            let r1 = unpack(&read[1].1[is-force_klen .. is]) ^ inverses[unpack(&read[0].1[0 .. force_klen])];
            let r2 = unpack(&read[1].2[is-force_klen .. is]) ^ !inverses[unpack(&read[0].2[0 .. force_klen])];
            let rscore = usize::count_ones(r0 & !r1 & !r2) as usize;   
            let r_ok = read[1].0[is .. is+klen].count_ones() == klen;

            let x0 = unpack(&read[0].0[is .. is + 9]) & unpack(&read[1].0[is .. is + 9]);
            let x1 = unpack(&read[0].1[is .. is + 9]) ^ unpack(&read[1].1[is .. is + 9]);
            let x2 = unpack(&read[0].2[is .. is + 9]) ^ unpack(&read[1].2[is .. is + 9]);
            let xscore = usize::count_ones(x0 & !x1 & !x2);

            if qscore == force_klen && rscore == force_klen && xscore == 9 && q_ok && r_ok {
                adapters.push((
                    false, // denotes r1 adapter seq candidate
                    unpack(&read[0].1[is .. is+klen]),
                    unpack(&read[0].2[is .. is+klen]),
                ));
                adapters.push((
                    true, // denotes r2 adapter seq candidate
                    unpack(&read[1].1[is .. is+klen]),
                    unpack(&read[1].2[is .. is+klen]),
                ));
            }   
            continue;
        }   

        if is >= klen && is + klen <= read[0].0.len() {
            let q0 = unpack(&read[0].0[is-klen .. is]) & inverses[unpack(&read[1].0[0 .. klen])];
            let q1 = unpack(&read[0].1[is-klen .. is]) ^ inverses[unpack(&read[1].1[0 .. klen])];
            let q2 = unpack(&read[0].2[is-klen .. is]) ^ !inverses[unpack(&read[1].2[0 .. klen])];
            let qscore = usize::count_ones(q0 & !q1 & !q2) as usize;   
            let q_ok = read[0].0[is .. is+klen].count_ones() == klen;
            if qscore == klen && q_ok {
                adapters.push((
                    false, // denotes r1 adapter seq candidate
                    unpack(&read[0].1[is .. is+klen]),
                    unpack(&read[0].2[is .. is+klen]),
                ));
            }
        }

        if is >= klen && is + klen <= read[1].0.len() {
            let r0 = unpack(&read[1].0[is-klen .. is]) & inverses[unpack(&read[0].0[0 .. klen])];
            let r1 = unpack(&read[1].1[is-klen .. is]) ^ inverses[unpack(&read[0].1[0 .. klen])];
            let r2 = unpack(&read[1].2[is-klen .. is]) ^ !inverses[unpack(&read[0].2[0 .. klen])];
            let rscore = usize::count_ones(r0 & !r1 & !r2) as usize;   
            let r_ok = read[1].0[is .. is+klen].count_ones() == klen;
            if rscore == klen && r_ok {
                adapters.push((
                    true, // denotes r2 adapter seq candidate
                    unpack(&read[1].1[is .. is+klen]),
                    unpack(&read[1].2[is .. is+klen]),
                ));
            }
        }
    }

    adapters
}

fn judge(exist: usize, diff: usize, klen: usize, maxdiff: usize) -> bool {

    if exist >= klen + maxdiff {
        return diff <= maxdiff;
    }

    if exist >= klen + maxdiff - 2 {
        return diff + 1 <= maxdiff;
    }

    if exist >= klen + maxdiff - 4 {
        return diff + 2 <= maxdiff;
    }

    if exist >= klen + maxdiff - 6 {
        return diff + 3 <= maxdiff;
    }

    if exist >= klen + maxdiff - 8 {
        return diff + 4 <= maxdiff;
    }

    if exist >= klen + maxdiff - 10 {
        return diff + 5 <= maxdiff;
    }

    if exist >= klen + maxdiff - 12 {
        return diff + 6 <= maxdiff;
    }

    return false;

}


fn quick_calc_trim_len(klen: usize, maxdiff: usize, read: &ReadData, altread: &ReadData, 
    adapter_bit1: usize, adapter_bit2: usize, inverses: &Vec<usize>) -> usize {
    // construct target kmer that we should see at boundary between insert & adapter
    // len is 2k bases, 0..k bases are inferred *end of insert* from start
    // k..2k bases are start of actual adapter
    let mut target0 = (inverses[unpack(&altread.0[.. klen])] << klen) 
        + ((1 << klen) - 1);
    let mut target1 = (inverses[unpack(&altread.1[.. klen])] << klen) 
        + adapter_bit1;

    // note: for valid revcomp, bit2 is always !-ed as well as inverted, so that A<->T and G<->C
    let mut target2 = (!inverses[unpack(&altread.2[.. klen])] << klen)
        + adapter_bit2;
    
    let l = read.0.len();

    // special case for is == -1
    let exist = unpack(&read.0[..klen-1]) & target0;
    let diff = exist & ((unpack(&read.1[..klen-1]) ^ target1 ) | (unpack(&read.2[..klen-1]) ^ target2));
    if usize::count_ones(exist) as usize >= klen - 2 && usize::count_ones(diff) as usize <= 2 {
        return 0;
    } 

    let mask = (1 << (klen + maxdiff)) - 1;

    for is in 0 .. l {
        // is => insert size
        let query0 = unpack(&read.0[max(is, klen) - klen .. min(is + klen, l)]);
        let query1 = unpack(&read.1[max(is, klen) - klen .. min(is + klen, l)]);
        let query2 = unpack(&read.2[max(is, klen) - klen .. min(is + klen, l)]);

        if is + klen > l {
            // getting near the end of the read so start pruning bases from end of target
            target0 >>= 1;
            target1 >>= 1;
            target2 >>= 1;
        }

        let mut exist = query0 & target0;
        let mut diff = exist & ((query1 ^ target1) | (query2 ^ target2));

        let mut iter = 0;
        while iter <= 2 * klen - diff {
            if judge(usize::count_ones(exist & mask) as usize, usize::count_ones(diff & mask) as usize, klen, maxdiff) {
                // this is where we want to trim to
                return is;
            }

            iter += 1;
            exist >>= 1;
            diff >>= 1;
        }
    }

    l
}

fn build_inverses(klen: usize) -> Vec<usize> {
    // inversion: 0x11010001 <-> 0x10001011 (flips bit order)
    // we pre-generate the inverse of every 'kmer' from 0 to 2^klen
    let bitmap = (0..klen).into_iter().map(|i| ((1 << i), ((2 * i as isize) + 1 - klen as isize))).collect::<Vec<_>>();

    (0 as usize .. (1 << klen)).into_par_iter().map(|k| {
        bitmap.iter().map(|(bit, offset)| 
            if offset > &0 {(k & *bit) >> *offset} else {(k & *bit) << -offset}
        ).sum::<usize>()
    }).collect::<Vec<_>>()
}


fn consensus(candidates: &Vec<(bool, usize, usize)>, klen: usize, side: bool, allow_fail: bool) -> (bool, usize, usize, usize) {
    // merges all candidates, picks consensus bases
    // panics if there isn't an obvious consensus at some base
    let mut per_base_votes = vec![[0; 4]; klen];
    for (s, a1, a2) in candidates.iter() {
        if s != &side {
            continue;
        }
        for j in 0..klen {
            let mut base = 0;
            if (a1 >> j) & 1 == 1 {
                base += 2;
            }
            if (a2 >> j) & 1 == 1 {
                base += 1;
            }
            per_base_votes[klen-j-1][base] += 1;
        }
    }
    let mut valid = true;
    
    let adapter = (0..klen).into_iter().map(|j| {
        let ranked_decision = per_base_votes[j].iter().enumerate()
            .sorted_by_key(|x| -x.1).collect::<Vec<_>>();
        if ranked_decision[1].1 + ranked_decision[2].1 / 2 + ranked_decision[3].1 / 3 > *ranked_decision[0].1 {
            println!("Warning: adapter may be unreliable at position {j} {ranked_decision:?} {}!", candidates.len());
            valid = false;
            if !allow_fail {
                std::process::exit(95);
            }

        }
        ranked_decision[0].0 << (2 * (klen - 1 - j))
    }).sum::<usize>();

    let adapter_bit1 = (0..klen).into_iter().map(|i| 
        (((adapter >> (1+(2*i))) & 1) as usize) << i
    ).sum::<usize>();
    let adapter_bit2 = (0..klen).into_iter().map(|i| 
        (((adapter >> (2*i)) & 1) as usize) << i
    ).sum::<usize>();

    (valid, adapter, adapter_bit1, adapter_bit2)
}


fn trim(inpaths: [String; 2], outpaths: [String; 2], klen: usize, maxdiff: usize) {
    let inverses = build_inverses(klen);
    let t0 = SystemTime::now();

    let readset = load_readset(inpaths[0].as_str(), inpaths[1].as_str());
    let t1 = SystemTime::now();

    let n_bases = readset.par_iter().map(|y| y[0].0.len() + y[1].0.len()).sum::<usize>();
    println!("Loaded {n_bases} bases {} reads from {inpaths:?} -- {:?}.\n", 
        readset.len(), t1.duration_since(t0).unwrap());

    let candidates = readset.par_iter().flat_map(|read| {
        get_adapter_candidates(read, &inverses, klen, maxdiff)
    }).collect::<Vec<_>>();
    
    let t2 = SystemTime::now();
    println!("Found {:?} strong adapter candidates -- {:?}.", 
        candidates.len(), t2.duration_since(t1).unwrap());
    
    let allow_fail = candidates.len() < 100;
    let (valid1, adapter1, a1_bit1, a1_bit2) = consensus(&candidates, klen, false, allow_fail);
    let (valid2, adapter2, a2_bit1, a2_bit2) = consensus(&candidates, klen, true, allow_fail);

    let trimlens = if !valid1 || !valid2 {
        println!("Couldn't identify adapter, however only {} adapter candidates, so skipping...", candidates.len());
        readset.par_iter().map(|read| {
            return (0, 0);
        }).collect::<Vec<_>>()
    } else {
        println!("Deduced adapter for trimming read 1: {}.", render(adapter1, klen));
        println!("Deduced adapter for trimming read 2: {}.", render(adapter2, klen));

        readset.par_iter().map(|read| {
            if read[0].0.len() < 32 || read[1].0.len() < 32 {
                return (0, 0);
            }

            let trim1 = quick_calc_trim_len(klen, maxdiff, &read[0], &read[1], a1_bit1, a1_bit2, &inverses);
            let trim2 = quick_calc_trim_len(klen, maxdiff, &read[1], &read[0], a2_bit1, a2_bit2, &inverses);

            if trim1 == read[0].0.len() && trim2 == read[1].0.len() {
                // special case when already pre trimmed to different lengths
                return (0, 0);
            }

            let trim = min(trim1, trim2);

            (trim, read[0].0.len() + read[1].0.len() - trim*2) 
        }).collect::<Vec<_>>()
    };

    let t3 = SystemTime::now();

    // hackily edit the original fastq trings to remove bases
    let outputs = (0..2).into_iter().map(|i| {
        readset.par_iter().zip(trimlens.par_iter()).filter_map(|(read, (trimlen, trim))| {
            if read[0].0.len() < 32 || read[1].0.len() < 32 || (trimlen > &0 && trimlen < &32) {
                return None;
            }
            
            if trim == &0 {
                return Some((read[i].0.len(), 0, format!("{}\n{}\n+\n{}\n", read[i].3, read[i].4, read[i].5)));
            }
            return Some((*trimlen, read[i].0.len() - trimlen, format!("{}\n{}\n+\n{}\n", read[i].3, 
                read[i].4[..*trimlen].to_string(), 
                read[i].5[..*trimlen].to_string()
            ))); 
        }).collect::<Vec<_>>()
    }).collect::<Vec<_>>();

    let pre_nr = readset.len();
    let pre_nb = readset.iter().map(|r| r[0].0.len() + r[1].0.len()).sum::<usize>();
    let post_nr = outputs[0].iter().zip(outputs[1].iter()).map(|(r1, r2)| if r1.1 == 0 && r2.1 == 0 {1} else {0}).sum::<usize>();
    let post_nb = outputs[0].iter().map(|r| r.0).sum::<usize>() + outputs[1].iter().map(|r| r.0).sum::<usize>() ;
    let trimmed_nr = pre_nr - post_nr;
    let trimmed_nb = pre_nb - post_nb;

    outputs.into_par_iter().enumerate().for_each(|(i, s)| {
        let mut tmp = GzEncoder::new(Vec::new(), Compression::fast());
        let _ = tmp.write_all(s.into_iter().map(|x| x.2).join("").as_bytes()).unwrap();
        let _ = fs::write(&outpaths[i], tmp.finish().unwrap());
    });

    let log_content = format!("Input:\t{} reads\t{} bases.\nTotal Removed:\t{} reads ({:.2}%)\t{} bases ({:.2}%)\nResult:\t{} reads ({:.2}%)\t{} bases ({:.2}%)\n",
        pre_nr, pre_nb,
        trimmed_nr, 100.0*trimmed_nr as f64 / pre_nr as f64, 
        trimmed_nb, 100.0*trimmed_nb as f64 / pre_nb as f64,
        post_nr, 100.0*post_nr as f64 / pre_nr as f64, 
        post_nb, 100.0*post_nb as f64 / pre_nb as f64,
    );

    let _ =fs::write("bbduk.log", log_content.clone());

    println!("\n{log_content}");

    let t4 = SystemTime::now();
    println!("Written output fastqs & bbduk.log -- {:?}.", 
        t4.duration_since(t3).unwrap(),
    );
}


#[derive(Parser, Debug)]
#[command(version, about, long_about = None)]
struct TrimArgs {
    #[arg(long, default_value_t = 27)] // max of ~31 (or performance suffers)
    klen: usize,
    #[arg(long, default_value_t = 8)] // how many error bases to tolerate
    maxdiff: usize,
    #[arg(long)]
    r1: String,
    #[arg(long)]
    r2: String,
    #[arg(long)]
    or1: String,
    #[arg(long)]
    or2: String,
}

fn main() {
    let args = TrimArgs::parse();

    if !Path::new(args.r1.as_str()).exists() || !Path::new(args.r2.as_str()).exists() {
        panic!("{} or {} does not exist.", args.r1, args.r2);
    }

    trim([args.r1, args.r2], [args.or1, args.or2], args.klen, args.maxdiff);
}
