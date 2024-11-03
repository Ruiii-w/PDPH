#ifndef PDPH_H
#define PDPH_H

#include <cstring>
#include <random>
#include <set>
#include <iostream>
#include <string>
#include <map>
#include <thread>
#include <mutex>
#include <functional>
#include <future>
#include <queue>
#include <vector>
#include "Hash.h"
#include "container.h"

template <typename T>
class SafeQueue {
private:
    std::queue<T> queue_;
    std::mutex mutex_;

public:
    SafeQueue() = default;
    SafeQueue(SafeQueue&& other) = default;
    ~SafeQueue() = default;

    bool empty() {
        std::lock_guard<std::mutex> lock(mutex_);
        return queue_.empty();
    }

    int size() {
        std::lock_guard<std::mutex> lock(mutex_);
        return queue_.size();
    }

    void enqueue(T&& task) {
        std::lock_guard<std::mutex> lock(mutex_);
        queue_.emplace(std::forward<T>(task));
    }

    bool dequeue(T& task) {
        std::lock_guard<std::mutex> lock(mutex_);
        if (queue_.empty()) return false;
        task = std::move(queue_.front());
        queue_.pop();
        return true;
    }
};


class ThreadPool {
private:
    class Worker {
    private:
        int id_;
        ThreadPool* pool_;

    public:
        Worker(ThreadPool* pool, int id) : pool_(pool), id_(id) {}

        void operator()() {
            std::function<void()> task;
            while (!pool_->shutdown_) {
                std::unique_lock<std::mutex> lock(pool_->conditional_mutex_);
                pool_->conditional_lock_.wait(lock, [&] { return pool_->shutdown_ || !pool_->task_queue_.empty(); });
                
                if (pool_->shutdown_ && pool_->task_queue_.empty()) return;
                
                if (pool_->task_queue_.dequeue(task)) {
                    task();
                }
            }
        }
    };

    bool shutdown_;
    SafeQueue<std::function<void()>> task_queue_;
    std::vector<std::thread> threads_;
    std::mutex conditional_mutex_;
    std::condition_variable conditional_lock_;

public:
    explicit ThreadPool(int thread_count = 4) : shutdown_(false), threads_(thread_count) {}

    ~ThreadPool() { shutdown(); }

    void init() {
        for (size_t i = 0; i < threads_.size(); ++i) {
            threads_[i] = std::thread(Worker(this, i));
        }
    }

    void shutdown() {
        shutdown_ = true;
        conditional_lock_.notify_all();
        for (auto& thread : threads_) {
            if (thread.joinable()) thread.join();
        }
    }

    template <typename F, typename... Args>
    auto submit(F&& func, Args&&... args) -> std::future<decltype(func(args...))> {
        auto task = std::make_shared<std::packaged_task<decltype(func(args...))()>>(
            std::bind(std::forward<F>(func), std::forward<Args>(args)...));
        std::function<void()> wrapper_func = [task]() { (*task)(); };
        task_queue_.enqueue(std::move(wrapper_func));
        conditional_lock_.notify_one();
        return task->get_future();
    }
};

template <class T>
class PDPH : Hash<T> {
public:
    uint64_t total_keys_;
    uint64_t actual_keys_;
    double s_factor;
    int container_count_;
    uint64_t extension_count_;

private:
    Container** containers_;

public:
    PDPH(uint64_t size, double _s_factor)
        : total_keys_(size), s_factor(_s_factor), actual_keys_(0),
          extension_count_(0) {
        initialize_containers(size / CONTAINER_INIT_SIZE + 1);
        initialize_hash_functions();
    }

    PDPH(uint64_t size, double _s_factor, size_t container_size)
        : total_keys_(size), s_factor(_s_factor), actual_keys_(0),
          extension_count_(0) {
        initialize_containers(size / container_size + 1, container_size);
        initialize_hash_functions();
    }

    ~PDPH() {
        for (int i = 0; i < container_count_; ++i) {
            delete containers_[i];
        }
        delete[] containers_;
    }

    int Insert(T key, ValueType value) { return perform_insert(key, value); }
    int Insert(T key, const uint32_t& fp, ValueType value) { return perform_insert(key, fp, value); }
    int Insert(T key, ValueType value, bool) { return perform_insert(key, value); }

    bool query(T key, uint32_t fp, uint64_t cell_id, uint16_t offset) { return retrieve(key, fp, cell_id, offset); }
    bool query(T key, uint32_t fp) { return retrieve(key, fp); }
    bool query(T key) { return retrieve(key); }
    double query_sequential(T key, uint32_t fp, uint64_t cell_id, uint16_t offset) { return retrieve_sequential(key, fp, cell_id, offset); }
    bool query_random(T key, uint32_t fp, uint64_t cell_id, uint16_t offset) { return retrieve_random(key, fp, cell_id, offset); }
    ValueType Get(T key) { return retrieve_value(key); }
    ValueType Get(T key, bool) { return retrieve_value(key); }

    void optimize_storage() {
        for (int i = 0; i < container_count_; ++i) {
            for (int j = 0; j < (1 << containers_[i]->global_depth); ++j) {
                (*(containers_[i]->segments[j]))->optimize();
            }
        }
    }

private:
    void initialize_containers(int count, size_t container_size = CONTAINER_INIT_SIZE) {
        container_count_ = count;
        containers_ = new Container*[container_count_];
        for (int i = 0; i < container_count_; ++i) {
            containers_[i] = new Container(container_size, s_factor);
        }
    }

    void initialize_hash_functions() {
        std::set<int> unique_seeds;
        for (int i = 0; i < MAXIMUM_LEVEL; ++i) {
            uint32_t seed = generate_unique_seed(unique_seeds);
            seed_of_level[i] = seed;
            unique_seeds.insert(seed);
        }
    }

    uint32_t generate_unique_seed(const std::set<int>& seeds) const {
        uint32_t seed = rand() % MAX_PRIME32;
        while (seeds.find(seed) != seeds.end()) {
            seed = rand() % MAX_PRIME32;
        }
        return seed;
    }

    int perform_insert(T key, const uint32_t& fingerprint, ValueType value) {
        int container_idx = fingerprint % container_count_;
        int segment_idx = key & ((1 << containers_[container_idx]->global_depth) - 1);
        int result = (*(containers_[container_idx]->segments[segment_idx]))->Insert(key, value);
        
        handle_segment_operations(container_idx, segment_idx, result);
        return result;
    }

    int perform_insert(T key, ValueType value) {
        uint32_t fingerprint = hash_(&key, sizeof(KeyType), seed_of_level[MAXIMUM_LEVEL - 1]);
        return perform_insert(key, fingerprint, value);
    }

    void handle_segment_operations(int container_idx, int segment_idx, int result) {
        auto* segment = (*(containers_[container_idx]->segments[segment_idx]));
        if (result == -1) {
            if (containers_[container_idx]->global_depth == segment->local_depth) {
                ++extension_count_;
                containers_[container_idx]->extend();
                containers_[container_idx]->split(segment_idx);
            } else {
                containers_[container_idx]->split(segment_idx);
            }
            delete segment;
        }
    }

    bool retrieve(T key, uint32_t fingerprint) {
        int container_idx = fingerprint % container_count_;
        return (*(containers_[container_idx]->segments[key & (1 << containers_[container_idx]->global_depth) - 1]))->Get(key);
    }

    bool retrieve(T key) {
        uint32_t fingerprint = hash_(&key, sizeof(KeyType), seed_of_level[MAXIMUM_LEVEL - 1]);
        return retrieve(key, fingerprint);
    }

    double retrieve_sequential(T key, uint32_t fingerprint, uint64_t cell_id, uint16_t offset) {
        int container_idx = fingerprint % container_count_;
        return (*(containers_[container_idx]->segments[key & (1 << containers_[container_idx]->global_depth) - 1]))->Get_SEQ(key);
    }

    bool retrieve_random(T key, uint32_t fingerprint, uint64_t cell_id, uint16_t offset) {
        int container_idx = fingerprint % container_count_;
        return (*(containers_[container_idx]->segments[key & (1 << containers_[container_idx]->global_depth) - 1]))->Get_RND(key);
    }

    ValueType retrieve_value(T key) {
        uint32_t fingerprint = hash_(&key, sizeof(KeyType), seed_of_level[MAXIMUM_LEVEL - 1]);
        return retrieve(key, fingerprint);
    }
};

#endif // PDPH_H
