#include <algorithm>
#include <fstream>
#include <gflags/gflags.h>
#include <string>
#include <sys/time.h>
#include <thread>
#include <unistd.h>
#include "./src/PDPH.h"

namespace param {
    std::string index_identifier, dataset_identifier, value_format;
    uint32_t thread_count;
}

DEFINE_string(index, "PDPHFilter", "Hash Index");
DEFINE_string(dataset, "ycsb_wl_0100_uniform.dat", "Dataset name");
DEFINE_string(fv, "fixed", "Fixed or variable key-value pairs");
DEFINE_uint32(thread, 1, "Number of threads");

void assign_cpu_affinity(uint32_t thread_index) {
    cpu_set_t cpu_set;
    CPU_ZERO(&cpu_set);
    CPU_SET((thread_index % 2 == 0) ? thread_index : thread_index + 1, &cpu_set);
    sched_setaffinity(0, sizeof(cpu_set_t), &cpu_set);
}

void concurrent_operation_on_extendible_ph(DataRange* thread_range, Hash<KeyType>* hash_table, std::vector<double>* latency_records = nullptr) {
    assign_cpu_affinity(thread_range->index);
    const uint64_t operation_start_index = thread_range->begin;
    const uint64_t operation_end_index = thread_range->end;
    const std::vector<DataRecord>& record_dataset = *thread_range->workload;

    int insertion_failed_count = 0, insertion_successful_count = 0, not_found_count = 0;
    int insertion_operation_count = 0, retrieval_operation_count = 0, update_operation_count = 0, default_operation_count = 0;

    struct timespec operation_start_time = {0, 0}, operation_end_time = {0, 0};

    for (uint64_t operation_index = operation_start_index; operation_index < operation_end_index; ++operation_index) {
        const DataRecord& current_record = record_dataset[operation_index];
        clock_gettime(CLOCK_REALTIME, &operation_start_time);
        
        const KeyType current_key = current_record.key;
        const uint32_t fingerprint = hash_(&current_key, sizeof(KeyType), seedOfLevel[MAXIMUM_LEVEL - 1]);
        const uint64_t cell_id = MurmurHash3_x86_32(&current_key, sizeof(KeyType), seed1);
        const uint16_t cell_offset = MurmurHash3_x86_32(&current_key, sizeof(KeyType), seed2) & cellMask;

        switch (current_record.op) {
            case DataRecord::INSERT_OPERATION:
                if (hash_table->Insert(current_record.key, current_record.value)) {
                    ++insertion_successful_count;
                } else {
                    ++insertion_failed_count;
                }
                ++insertion_operation_count;
                break;

            case DataRecord::GET_OPERATION:
                if (!hash_table->query(current_record.key, fingerprint, cell_id, cell_offset)) {
                    ++not_found_count;
                }
                ++retrieval_operation_count;
                break;

            default:
                ++default_operation_count;
        }

        clock_gettime(CLOCK_REALTIME, &operation_end_time);
        const double operation_latency = (double)(operation_end_time.tv_nsec - operation_start_time.tv_nsec) +
                                         (double)(operation_end_time.tv_sec - operation_start_time.tv_sec) * 1e9;
        if (latency_records) {
            latency_records->emplace_back(operation_latency);
        }
    }

    gettimeofday(&thread_range->tv, nullptr);
    std::cout << "Insertion failed: " << insertion_failed_count << " | Insertion count: " << insertion_operation_count << std::endl;
    std::cout << "Search not found: " << not_found_count << " | Search count: " << retrieval_operation_count << std::endl;
}

void benchmark_trace_on_extendible_ph(std::vector<DataRecord>* record_data, const size_t& total_records, const uint32_t& thread_count, void (*test_function)(DataRange*, Hash<KeyType>*, std::vector<double>*)) {
    const uint32_t trace_count = total_records / thread_count;
    std::vector<std::thread> thread_array(thread_count);
    std::vector<DataRange> thread_ranges(thread_count);
    timeval start_time, end_time;
    double total_duration = 0;

    Hash<KeyType>* hash_table = new PDPH<KeyType>(total_records, DEFAULT_SLT);

    for (uint32_t i = 0; i < thread_count; ++i) {
        thread_ranges[i] = {i, i * trace_count, (i + 1) * trace_count, rand(), 8, record_data};
        thread_array[i] = std::thread(test_function, &thread_ranges[i], hash_table);
    }

    gettimeofday(&start_time, nullptr);
    for (auto& thread : thread_array) {
        thread.join();
    }

    double max_duration = (double)(thread_ranges[0].tv.tv_usec - start_time.tv_usec) / 1e6 +
                          (double)(thread_ranges[0].tv.tv_sec - start_time.tv_sec);
    double min_duration = max_duration;
    total_duration = max_duration;

    for (uint32_t i = 1; i < thread_count; ++i) {
        double interval = (double)(thread_ranges[i].tv.tv_usec - start_time.tv_usec) / 1e6 +
                         (double)(thread_ranges[i].tv.tv_sec - start_time.tv_sec);
        total_duration += interval;
        min_duration = std::min(min_duration, interval);
        max_duration = std::max(max_duration, interval);
    }

    total_duration /= thread_count;
    std::cout << thread_count << " Threads, Total Time = " << total_duration << " s, Throughput = " << (record_data->size()) / total_duration / 1e6
              << " Mops/s" << std::endl;
}

void load_dataset_from_file(const std::string& file_path, std::vector<DataRecord>* record_data) {
    std::cout << "Loading Dataset from: " << file_path << std::endl;
    std::ifstream file_stream{file_path, std::ios::binary | std::ios::ate};
    if (!file_stream) {
        std::cerr << "Cannot open file: " << file_path << std::endl;
        return;
    }

    const auto file_size = file_stream.tellg();
    if (file_size == 0) {
        std::cerr << "Empty file" << std::endl;
        return;
    }

    const uint64_t record_count = file_size / sizeof(DataRecord);
    record_data->resize(record_count);
    
    if (!file_stream.read(reinterpret_cast<char*>(record_data->data()), file_size)) {
        std::cerr << "Error reading from " << file_path << std::endl;
    } else {
        std::cout << "Loading completed. " << record_data->size() << " key-value pairs loaded." << std::endl;
    }
}

void initialize_index_on_extendible_ph(const std::string& index_name, const std::string& dataset_name, const std::string& value_format, uint64_t thread_count) {
    std::vector<DataRecord> record_dataset;
    generateRandomKeyValueRecords(&record_dataset, 140000000);
    const uint64_t total_operations = record_dataset.size();
    std::cout << "============================= " << index_name << " =============================" << std::endl;

    if (index_name == "PDPHFilter") {
        benchmark_trace_on_extendible_ph(&record_dataset, total_operations, thread_count, concurrent_operation_on_extendible_ph);
    }
}

int main(int argc, char* argv[]) {
    google::ParseCommandLineFlags(&argc, &argv, true);
    param::dataset_identifier = FLAGS_dataset;
    param::thread_count = FLAGS_thread;
    param::value_format = FLAGS_fv;
    param::index_identifier = FLAGS_index;

    initialize_index_on_extendible_ph(param::index_identifier, param::dataset_identifier, param::value_format, param::thread_count);
}
