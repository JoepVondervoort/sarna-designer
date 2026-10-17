#include <array>
#include <cstdint>
#include <cstdio>
#include <iostream>
#include <string>
#include <vector>

// Independent rolling 2-bit exhaustive Hamming counter. No Bowtie output, seeding,
// index, annotation, or Python screening implementation is used to determine matches.
int main(int argc, char** argv) {
    if (argc < 3) return 2;
    constexpr int length = 19;
    constexpr uint64_t mask = (uint64_t(1) << (2 * length)) - 1;
    constexpr uint64_t lowbits = 0x5555555555555555ULL & mask;
    auto base = [](unsigned char c) -> int {
        switch (c) {case 'A': case 'a': return 0; case 'C': case 'c': return 1;
                    case 'G': case 'g': return 2; case 'T': case 't': return 3; default: return -1;}
    };
    std::vector<uint64_t> queries;
    for (int i = 2; i < argc; ++i) {
        std::string sequence(argv[i]);
        if (sequence.size() != length) return 3;
        uint64_t packed = 0;
        for (unsigned char c : sequence) {
            int value = base(c);
            if (value < 0) return 4;
            packed = (packed << 2) | uint64_t(value);
        }
        queries.push_back(packed);
    }
    std::vector<std::array<uint64_t,4>> counts(queries.size());
    for (auto& count : counts) count.fill(0);
    FILE* input = std::fopen(argv[1], "rb");
    if (!input) return 5;
    std::vector<unsigned char> buffer(4 * 1024 * 1024);
    uint64_t rolling = 0;
    int valid = 0;
    bool header = false, header_end = false;
    std::string seqid;
    uint64_t position = 0;
    unsigned exact_printed = 0;
    size_t received;
    while ((received = std::fread(buffer.data(), 1, buffer.size(), input)) != 0) {
        for (size_t i = 0; i < received; ++i) {
            unsigned char c = buffer[i];
            if (c == '>') { header = true; header_end = false; seqid.clear(); position=0; valid = 0; rolling = 0; continue; }
            if (c == '\n') {header = false; continue;}
            if (header) {
                if (c==' ' || c=='\t' || c=='\r') header_end=true;
                if (!header_end) seqid += c;
                continue;
            }
            if (c=='\r' || c==' ' || c=='\t') continue;
            ++position;
            int value = base(c);
            if (value < 0) {valid = 0; rolling = 0; continue;}
            rolling = ((rolling << 2) | uint64_t(value)) & mask;
            if (valid < length) ++valid;
            if (valid < length) continue;
            for (size_t q = 0; q < queries.size(); ++q) {
                uint64_t difference = rolling ^ queries[q];
                unsigned mismatches = __builtin_popcountll((difference | (difference >> 1)) & lowbits);
                if (mismatches <= 3) ++counts[q][mismatches];
                if (mismatches==0 && q>=4 && exact_printed++<100)
                    std::cerr << q << "," << seqid << "," << position-length << "\n";
            }
        }
    }
    if (std::ferror(input)) {std::fclose(input); return 6;}
    std::fclose(input);
    for (size_t q = 0; q < queries.size(); ++q)
        std::cout << q << ',' << counts[q][0] << ',' << counts[q][1] << ','
                  << counts[q][2] << ',' << counts[q][3] << '\n';
}
