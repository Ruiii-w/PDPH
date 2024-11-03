#ifndef UTILITIES_HEADER_FILE
#define UTILITIES_HEADER_FILE

#include <atomic>
#include <mutex>
#include <sys/stat.h>
#include <map>
#include <string>
#include <cmath>
#include "murmur3.h"
#include <vector>
#include <cstdint>
#include <cstdlib>
#include <immintrin.h>

#define LIKELY_CONDITION(x)       __builtin_expect((x),1)
#define UNLIKELY_CONDITION(x)     __builtin_expect((x),0)
#define KEY_DATA_TYPE uint64_t
#define VALUE_DATA_TYPE uint8_t
#define CONTAINER_SIZE N
#define CONTAINER_INIT_SIZE (1 << 16)
#define CONTAINER_OFF_MASK  (CONTAINER_INIT_SIZE - 1)
#define CONTAINER_POOL_SIZE 100
#define IS_BIT_SET(var, pos) ((((var) & (1 << (pos))) > 0) ? (1) : (0))
#define FREE_ARRAY_SPACE(arrayPtr) {if(arrayPtr != nullptr) delete[] arrayPtr; arrayPtr=nullptr;}
#define HASH_FUNCTION(key, length, seed) MurmurHash3_x86_32(key, length, seed)
#define MAXIMUM_LEVEL 30
#define CURRENT_LEVEL MAXIMUM_LEVEL
#define CAS_OPERATION(pointer, expected, desired) \
  (__atomic_compare_exchange_n(pointer, expected, desired, false, __ATOMIC_ACQUIRE, \
                                __ATOMIC_ACQUIRE))

using KeyType = size_t;
using ValueType = VALUE_DATA_TYPE;
const uint32_t lock_pos = ((uint32_t)1 << 31);
const uint32_t lock_mask = ((uint32_t)1 << 31) - 1;
const uint32_t fp_mask = (1 << 8) - 1;
const uint64_t key_mask = ((uint64_t)1 << 8)  - 1;
const KeyType SENTINEL_KEY = static_cast<KeyType>(-2);  
const KeyType INVALID_KEY = static_cast<KeyType>(-1);   
const ValueType NONE_VALUE = static_cast<ValueType>(0x0);
const uint64_t seed1 = 101, seed2 = 102, seed3 = 103, seed4=200;
const uint64_t cellSize[MAXIMUM_LEVEL] = {8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8};
const uint64_t offsetLength[MAXIMUM_LEVEL] = {6, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5};
const double key_to_cell_rate = 0.1;
const uint64_t cellnum0mask = ((1<<24)-1), totalSlotMask=(1ull<<27)-1;
const uint64_t cellpipe = (1 << 14), cellMask = cellpipe - 1;

struct DynamicStringKey {
    int length;
    char key[0];
};

template <typename MutexType>
class ScopedLock {
public:
    explicit ScopedLock(MutexType& mutex) : m_mutex(mutex) {
        m_mutex.lock();
    }

    ~ScopedLock() {
        m_mutex.unlock();
    }

private:
    MutexType &m_mutex;
};

template<typename KeyType>
struct KeyValuePair {
    KeyType key;
    ValueType value;

    KeyValuePair() : key{INVALID_KEY} {}
    KeyValuePair(KeyType _key, ValueType _value) : key{_key}, value{_value} {}

    KeyValuePair& operator=(const KeyValuePair& other) {
        key = other.key;
        value = other.value;
        return *this;
    }

    void* operator new(size_t size) {
        void* ptr;
        posix_memalign(&ptr, 64, size);
        return ptr;
    }

    void* operator new[](size_t size) {
        void* ptr;
        posix_memalign(&ptr, 64, size);
        return ptr;
    }
};


uint32_t seed_of_level[MAXIMUM_LEVEL];
enum class ContainerCategory {
    BucketMapContainer, 
    PerfectHashContainer
};

struct DataRecord {
    enum OperationType : uint32_t {INSERT_OPERATION = 0, GET_OPERATION = 1, UPDATE_OPERATION = 2};
    OperationType operation = INSERT_OPERATION;
    uint64_t key;       
    uint64_t value;
    uint64_t cell_id;
    uint32_t fingerprint;
    uint16_t offset;

    bool operator==(const LogRecord& other) const {
        return key == other.key && operation == other.operation;
    }

    bool operator<(const LogRecord& other) const {
        return key < other.key;
    }
};

struct DataRange {
    int index;
    uint64_t begin;
    uint64_t end;
    int length;
    std::vector<DataRecord> *workload;
    uint64_t random_number;
    struct timeval timestamp;
};

static bool IsFileExists(const char *file_path) {
    struct stat buffer;
    return (stat(file_path, &buffer) == 0);
}

template<typename ValueType>
class MemoryCell {
private:
    std::atomic<uint32_t> version_lock; 
    std::atomic<uint8_t> local_depth;   
    ValueType value;                    
    std::mutex mtx;                

public:
    MemoryCell() : version_lock(0), local_depth(0), value{} {}

    
    ValueType GetValue() {
        std::lock_guard<std::mutex> guard(mtx);
        return value;
    }

    
    void SetValue(const ValueType& new_value) {
        std::lock_guard<std::mutex> guard(mtx);
        value = new_value;
        version_lock.fetch_add(1, std::memory_order_relaxed); 
    }

    
    bool TryLock() {
        uint32_t expected = version_lock.load(std::memory_order_acquire);
        return version_lock.compare_exchange_strong(expected, expected + 1, std::memory_order_acq_rel);
    }

    
    void Unlock() {
        version_lock.fetch_sub(1, std::memory_order_release);
    }

    
    bool ValidateLockSet(uint32_t &expected_version) {
        return version_lock.load(std::memory_order_acquire) == expected_version;
    }

    
    void SetLocalDepth(uint8_t depth) {
        std::lock_guard<std::mutex> guard(mtx);
        local_depth.store(depth, std::memory_order_relaxed);
    }

    uint8_t GetLocalDepth() {
        return local_depth.load(std::memory_order_relaxed);
    }
};

static bool file_exists(const char *pool_path)
{
    struct stat buffer;
    return (stat(pool_path, &buffer) == 0);
}

void generateRandomKeyValueRecords(std::vector<Record>* keyValueStore, int recordCount) 
{
    for (int recordIndex = 1; recordIndex <= recordCount; ++recordIndex) 
    {
        Record recordEntry;
        recordEntry.key = MurmurHash3_x86_32(&recordIndex, sizeof(int), hashSeedLevel3);
        recordEntry.value = recordIndex;
        recordEntry.operationType = Record::Operation::INSERT;

        keyValueStore->emplace_back(recordEntry);
    }
}


#endif // UTILITIES_HEADER_FILE