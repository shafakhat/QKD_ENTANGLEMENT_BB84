# Drop the counts you downloaded from Classiq Studio here

Name the files after the circuit they came from, so the analyser can find them:

    chsh_a0b0.json   chsh_a0b1.json   chsh_a1b0.json   chsh_a1b1.json
    qrng_rotation.json
    qrng_ghz.json
    bb84_batch.json

.json or .csv, any of these layouts:

    {"00": 501, "11": 499}                      plain dict
    {"counts": {"00": 501, "11": 499}}          wrapped
    [{"bitstring": "00", "counts": 501}, ...]   list of records
    bitstring,counts                            CSV columns

Then:

    python3 code/analyze_studio_counts.py

The four CHSH files give you a Classiq-measured S value. `qrng_rotation` gives
P(1), the bias and the min-entropy of the Classiq samples.
