use std::env;
use std::fs::File;
use std::io::Read;
use bio::io::fasta;
use bio::scores::blosum62;
use bio::alignment::{pairwise, AlignmentOperation};
use phylogeny::{neighbor_joining, DistanceMatrix};

fn main() {
    // Parse command line arguments
    let args: Vec<String> = env::args().collect();
    if args.len() != 3 {
        panic!("Please provide a sequence as argument");
    }
    let input_file = &args[1];
    let target_database = &args[2];

    // Read the input sequence
    let mut file = File::open(input_file).expect("Unable to open file");
    let mut contents = String::new();
    file.read_to_string(&mut contents).expect("Unable to read file");
    let reader = fasta::Reader::new(contents.as_bytes());
    let input_record = reader.records().next().expect("Unable to read record").expect("Unable to unwrap record");
    let input_sequence = input_record.seq().to_owned();

    // Read the "Viola_philippica.aa" file
    let mut file = File::open(target_database).expect("Unable to open file");
    let mut contents = String::new();
    file.read_to_string(&mut contents).expect("Unable to read file");
    let reader = fasta::Reader::new(contents.as_bytes());

    // Align the sequences
    let gap_open = -5;
    let gap_extend = -1;

    let mut sequences = Vec::new();
    for record in reader.records() {
        let record = record.expect("Unable to read record");
        // get the sequence and id from record
        let viola_sequence = record.seq().to_owned();
        let id = record.id().to_owned();

        let mut aligner = pairwise::Aligner::with_capacity(input_sequence.len(), viola_sequence.len(), gap_open, gap_extend, &blosum62);

        let alignment = aligner.local(&input_sequence, &viola_sequence);

        // Calculate sequence identity by counting the number of matching positions and put its sequence name and sequence to sequences vector
        let numer_of_matches = alignment.operations.iter()
        .filter(|&op| *op == AlignmentOperation::Match).count();
        let percent_identity = numer_of_matches as f64 / alignment.operations.len() as f64;
        if percent_identity > 0.8 { // TODO: change 0.8 to a command line argument
            sequences.push((id, viola_sequence));
        }
    }

    // Create a distance matrix
    let mut distance = Vec::new();
    for (id, sequence) in &sequences {
        let mut row = Vec::new();
        for (id2, sequence2) in &sequences {
            let mut aligner = pairwise::Aligner::with_capacity(sequence.len(), sequence2.len(), gap_open, gap_extend, &blosum62);
            let alignment = aligner.local(&sequence, &sequence2);
            let numer_of_matches = alignment.operations.iter()
            .filter(|&op| *op == AlignmentOperation::Match).count();
            let percent_identity = numer_of_matches as f32 / alignment.operations.len() as f32;
            row.push(1.0 - percent_identity);
        }
        distance.push(row);
    }

    let distance_matrix = DistanceMatrix::new(sequences.iter().map(|(id, _)| id.clone()).collect(), distance);

    let tree = neighbor_joining(&distance_matrix);
    

}