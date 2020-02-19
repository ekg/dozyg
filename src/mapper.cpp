#include "mapper.hpp"

namespace gyeet {

using namespace std::chrono_literals;

void for_each_seq_in_file(
    const std::string& filename,
    const std::function<void(const std::string&, const std::string&)>& func) {
    // detect file type
    bool input_is_fasta = false;
    bool input_is_fastq = false;
    std::string line;
    igzstream in(filename.c_str());
    std::getline(in, line);
    if (line[0] == '>') {
        input_is_fasta = true;
    } else if (line[0] == '@') {
        input_is_fastq = true;
    } else {
        std::cerr << "unknown file format given to seqindex_t" << std::endl;
        assert(false);
        exit(1);
    }
    if (input_is_fasta) {
        while (in.good()) {
            std::string name = line.substr(1, line.find(" "));
            std::string seq;
            while (std::getline(in, line)) {
                if (line[0] == '>') {
                    // this is the header of the next sequence
                    break;
                } else {
                    seq.append(line);
                }
            }
            func(name, seq);
        }
    } else if (input_is_fastq) {
        while (in.good()) {
            std::string name = line.substr(1, line.find(" "));
            std::string seq;
            std::getline(in, seq); // sequence
            std::getline(in, line); // delimiter
            std::getline(in, line); // quality
            std::getline(in, line); // next header
            func(name, seq);
        }
    }
}
                          
void reader_thread(
    const std::vector<std::string>& files,
    seq_atomic_queue_t& seq_queue,
    std::atomic<bool>& reader_done) {
    for (auto& file : files) {
        for_each_seq_in_file(
            file,
            [&seq_queue](
                const std::string& name,
                const std::string& seq) {
                seq_record_t* rec = new seq_record_t(name, seq);
                while (!seq_queue.try_push(rec)) {
                    std::this_thread::sleep_for(100ns);
                }
            });
        
    }
    reader_done.store(true);
}

void writer_thread(
    std::ostream& out,
    gaf_atomic_queue_t& gaf_queue,
    const std::vector<std::atomic<bool>>& working) {
    while (true) {
        std::string* gaf_lines = nullptr;
        if (!gaf_queue.try_pop(gaf_lines)
            && !still_working(working)) {
            break;
        } else if (gaf_lines != nullptr) {
            out << *gaf_lines;
            delete gaf_lines;
        } else {
            std::this_thread::sleep_for(100ns);
        }
    }
}

bool still_working(
    const std::vector<std::atomic<bool>>& working) {
    bool ongoing = false;
    for (auto& w : working) {
        ongoing = ongoing || w.load();
    }
    return ongoing;
}

std::string map_seq(
    const std::string& name,
    const std::string& query,
    const gyeet_index_t& index,
    const uint64_t& max_chain_gap,
    const double& mismatch_rate,
    const uint64_t& align_best_n,
    const bool& write_alignments,
    const bool& write_chains,
    const bool& write_superchains) {

    auto anchors = anchors_for_query(index,
                                     query.c_str(),
                                     query.length());
    auto query_chains = chains(anchors,
                               index.kmer_length,
                               max_chain_gap,
                               mismatch_rate);
    std::stringstream ss;
    if (write_chains) {
        for (auto& chain : query_chains) {
            write_chain_gaf(ss, chain, index, name, query.length());
        }
    }
    auto query_superchains = superchains(query_chains, index.kmer_length, mismatch_rate);
    if (write_superchains) {
        for (auto& superchain : query_superchains) {
            write_superchain_gaf(ss, superchain, index, name, query.length());
        }
    }
    uint64_t up_to = std::min(align_best_n, (uint64_t)query_superchains.size());
    if (write_alignments) {
        //std::cerr << "aligning " << name << std::endl;
        for (uint64_t i = 0; i < up_to; ++i) {
            auto& superchain = query_superchains[i];
            alignment_t aln = superalign(
                name,
                query.length(),
                query.c_str(),
                superchain,
                index,
                index.kmer_length,
                mismatch_rate,
                max_chain_gap);
            write_alignment_gaf(ss, aln, index);
        }
    }
    return ss.str();
}

void worker_thread(
    uint64_t tid,
    std::atomic<bool>& is_working,
    seq_atomic_queue_t& seq_queue,
    gaf_atomic_queue_t& gaf_queue,
    const std::atomic<bool>& reader_done,
    const gyeet_index_t& index,
    uint64_t max_chain_gap,
    double mismatch_rate,
    uint64_t align_best_n,
    bool write_alignment,
    bool write_chains,
    bool write_superchains) {

    is_working.store(true);
    while (true) {
        seq_record_t* rec = nullptr;
        if (!seq_queue.try_pop(rec)
            && reader_done.load()) {
            break;
        } else if (rec != nullptr) {
            std::string* gaf_rec
                = new std::string(
                    map_seq(rec->name,
                            rec->seq,
                            index,
                            max_chain_gap,
                            mismatch_rate,
                            align_best_n,
                            write_alignment,
                            write_chains,
                            write_superchains));
            gaf_queue.push(gaf_rec);
            delete rec;
        } else {
            std::this_thread::sleep_for(100ns);
        }
    }
    is_working.store(false);
}

void map_reads(
    const std::vector<std::string>& input_files,
    const gyeet_index_t& index,
    const uint64_t& max_chain_gap,
    const double& mismatch_rate,
    const uint64_t& align_best_n,
    const uint64_t& nthreads,
    const bool& write_alignment,
    const bool& write_chains,
    const bool& write_superchains) {

    seq_atomic_queue_t seq_queue;
    gaf_atomic_queue_t gaf_queue;
    std::atomic<bool> reader_done;
    reader_done.store(false);
    std::vector<std::atomic<bool>> working(nthreads);
    for (auto& w : working) {
        w.store(true);
    }
    // launch reader
    std::thread reader(&reader_thread, std::ref(input_files), std::ref(seq_queue), std::ref(reader_done));
    // launch writer
    std::thread writer(&writer_thread, std::ref(std::cout), std::ref(gaf_queue), std::ref(working));
    // launch workers
    std::vector<std::thread> workers; workers.reserve(nthreads);
    for (uint64_t t = 0; t < nthreads; ++t) {
        workers.emplace_back(&worker_thread,
                             t,
                             std::ref(working[t]),
                             std::ref(seq_queue),
                             std::ref(gaf_queue),
                             std::ref(reader_done),
                             std::ref(index),
                             max_chain_gap,
                             mismatch_rate,
                             align_best_n,
                             write_alignment,
                             write_chains,
                             write_superchains);
    }

    reader.join();
    //std::cerr << "reader done" << std::endl;
    for (auto& worker : workers) {
        worker.join();
    }
    //std::cerr << "workers done" << std::endl;
    writer.join();
    //std::cerr << "writer done" << std::endl;
    //while (!reader_done.load() || !still_working(working)) {
//}
    // watch until we're done reading, done mapping, and done writing
    // then join and finish
    
}

}
