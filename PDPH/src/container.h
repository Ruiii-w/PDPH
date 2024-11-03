#ifndef CONTAINER_H
#define CONTAINER_H
#include "util.h"
#include <atomic>
#include <cstring>
#include <memory>
#include <mutex>
#include <iostream>
#include <sys/time.h>
#include <unordered_map>

#define LOCK_CLOCK
#define DEFAULT_SLT 0.85
#define DEFAULT_TRANS 0.5
#define DEFAULT_MF 10

const uint8_t BIT_SHIFT_MAST = 0x0f;
double part0 = 0, part1 = 0, part2 = 0;
struct timeval t1,t2,t3,t4,t5,t6;


class BitArray {
private:
    uint32_t multiple_factor;
    uint32_t bitArraySize;
    uint8_t* bitArray;

public:
    BitArray(uint32_t _multiple_factor = DEFAULT_MF) 
        : multiple_factor(_multiple_factor), bitArraySize(0), bitArray(nullptr) {
        initializeBitArray();
    }

    void initializeBitArray() {
        bitArraySize = ((CONTAINER_INIT_SIZE * multiple_factor) >> 3) + 1;
        bitArray = new uint8_t[bitArraySize]();
    }

    void setBit(uint32_t index) {
        uint32_t byteIndex = index >> 3;
        uint32_t bitIndex = index & 0x7;
        if (byteIndex >= bitArraySize) {
            byteIndex %= bitArraySize; 
        }
        bitArray[byteIndex] |= (1 << bitIndex);
    }

    bool isBitSet(uint32_t index) const {
        uint32_t byteIndex = index >> 3;
        uint32_t bitIndex = index & 0x7;
        if (byteIndex >= bitArraySize) {
            byteIndex %= bitArraySize; 
        }
        return (bitArray[byteIndex] & (1 << bitIndex)) != 0;
    }

    void clearBit(uint32_t index) {
        uint32_t byteIndex = index >> 3;
        uint32_t bitIndex = index & 0x7;
        if (byteIndex >= bitArraySize) {
            byteIndex %= bitArraySize; 
        }
        bitArray[byteIndex] &= ~(1 << bitIndex);
    }

    void clearAllBits() {
        std::fill_n(bitArray, bitArraySize, 0);
    }

    uint32_t getSize() const {
        return bitArraySize;
    }

    void dumpToFile(const char* filePath) const {
        FILE* filePointer = fopen(filePath, "wb");
        if (filePointer) {
            fwrite(bitArray, sizeof(uint8_t), bitArraySize, filePointer);
            fclose(filePointer);
        }
    }

    void loadFromFile(const char* filePath) {
        FILE* filePointer = fopen(filePath, "rb");
        if (filePointer) {
            fread(bitArray, sizeof(uint8_t), bitArraySize, filePointer);
            fclose(filePointer);
        }
    }

    uint32_t extract32BitSlice(uint32_t offset) const {
        uint32_t ret = 0;
        uint32_t byteIndex = offset >> 3;
        uint32_t bitIndex = offset & 0x7;
        uint32_t sliceBuffer[5] = {0};

        for (int i = 0; i < 5; ++i) {
            uint32_t currentByteIndex = byteIndex + i;
            if (currentByteIndex >= bitArraySize) {
                currentByteIndex %= bitArraySize; 
            }
            sliceBuffer[i] = bitArray[currentByteIndex];
        }

        ret = (sliceBuffer[4] << (32 - bitIndex)) | 
              (sliceBuffer[3] << (24 - bitIndex)) | 
              (sliceBuffer[2] << (16 - bitIndex)) | 
              (sliceBuffer[1] << (8 - bitIndex)) | 
              (sliceBuffer[0] >> bitIndex);
        
        return ret;
    }

    ~BitArray() {
        delete[] bitArray;
    }
};


class TransHash {
private:
    int effectiveLevelLimit;
    uint64_t totalKeys, totalSlots, totalCells;
    uint64_t cellCounts[MAXIMUM_LEVEL];
    uint8_t* cellArray[MAXIMUM_LEVEL];
    std::unordered_map<KeyType, int> hierarchyLevels;
    uint64_t* cellKeyCounts;
    uint32_t* bitVector;
    KeyType* primaryKeyArray, *primaryKeyCache;
    ValueType* primaryValueArray, *primaryValueCache;
    uint64_t memoryBitUsage;

public:
    explicit TransHash() {
        memset(cellCounts, 0, sizeof(cellCounts));
        memset(cellArray, 0, sizeof(cellArray));
    }

    bool initialize(KeyType* keys, ValueType* primaryValueCache, uint64_t keyCount, uint64_t slotCount, uint64_t cellCounts) {
        memoryBitUsage = 0;
        effectiveLevelLimit = 0;
        totalKeys = keyCount;
        totalSlots = ((slotCount >> 5) + 1) << 5;
        totalCells = cellCounts;

        cellKeyCounts = new uint64_t[totalCells + 1];
        bitVector = new uint32_t[totalSlots >> 5];
        memset(bitVector, 0, (totalSlots >> 5) * sizeof(uint32_t));

        primaryKeyArray = new KeyType[totalSlots];
        memset(primaryKeyArray, 0, totalSlots * sizeof(KeyType));
        primaryValueArray = new ValueType[totalSlots];
        memset(primaryValueArray, 0, totalSlots * sizeof(ValueType));
        primaryKeyCache = new KeyType[totalSlots];
        memset(primaryKeyCache, 0, totalSlots * sizeof(KeyType));
        primaryValueCache = new ValueType[totalSlots];
        memset(primaryValueCache, 0, totalSlots * sizeof(ValueType));

        memoryBitUsage += totalSlots;

        initializeCell(keys, primaryValueCache, totalKeys, 0);
        return true;
    }

    void initialize(uint64_t keyCount, uint64_t slotCount, uint64_t cellCounts, uint32_t* bitMaskArray) {
        memoryBitUsage = 0;
        totalKeys = keyCount;
        totalSlots = (slotCount >> 5) << 5;
        totalCells = cellCounts;

        cellKeyCounts = new uint64_t[totalCells + 1];
        bitVector = bitMaskArray;
        primaryKeyArray = new KeyType[totalSlots];
        memset(primaryKeyArray, 0, totalSlots * sizeof(KeyType));
        primaryValueArray = new ValueType[totalSlots];
        memset(primaryValueArray, 0, totalSlots * sizeof(ValueType));
        primaryKeyCache = new KeyType[totalSlots];
        memset(primaryKeyCache, 0, totalSlots * sizeof(KeyType));
        primaryValueCache = new ValueType[totalSlots];
        memset(primaryValueCache, 0, totalSlots * sizeof(ValueType));
    }

    void initializeCell(KeyType* keys, ValueType* primaryValueCache, uint64_t keyCount, int level) {
        if (level == 0) cellCounts[level] = keyCount * 0.17 + 1;
        if (level == 1) cellCounts[level] = keyCount * 0.46 + 1;
        if (level >= 2) cellCounts[level] = keyCount + 1001;

        memoryBitUsage += cellCounts[level] * cellSize[level];
        cellArray[level] = new uint8_t[cellCounts[level]];
        effectiveLevelLimit = std::max(effectiveLevelLimit, level);
        memset(cellKeyCounts, 0, (cellCounts[level] + 1) * sizeof(uint64_t));

        for (uint64_t i = 0; i < keyCount; i++) {
            KeyType key = keys[i];
            uint64_t cellID = MurmurHash3_x86_32(&key, sizeof(KeyType), seed1) % cellCounts[level];
            cellKeyCounts[cellID]++;
        }

        for (uint64_t j = 0, maxCount = cellCounts[level], lastCellCount = 0; j < maxCount; j++) {
            if (cellKeyCounts[j] == 0) {
                cellKeyCounts[j] = -1;
            } else {
                cellKeyCounts[j] += lastCellCount;
                lastCellCount = cellKeyCounts[j];
            }
        }
        for (uint64_t i = 0; i < keyCount; i++) {
            KeyType key = keys[i];
            uint64_t cellID = MurmurHash3_x86_32(&key, sizeof(KeyType), seed1) % cellCounts[level];
            uint64_t keyPosition = --cellKeyCounts[cellID];
            primaryKeyCache[keyPosition] = key;
            primaryValueCache[keyPosition] = primaryValueCache[i];
        }

        std::vector<KeyType> failedKeys;
        std::vector<ValueType> failedValues;
        uint64_t rangeStart[10] = {7, 6, 5, 4, 3, 2, 1};
        uint64_t rangeEnd[10] = {100, 6, 5, 4, 3, 2, 1};

        cellKeyCounts[cellCounts[level]] = keyCount;

        for (int t = 0; t < 7; t++) {
            uint64_t failedKeyCount = 0;
            for (uint64_t j = 0, maxCount = cellCounts[level]; j < maxCount; j++) {
                if (cellKeyCounts[j] == -1) continue;
                uint64_t k = j + 1;
                while (k <= maxCount && cellKeyCounts[k] == -1) k++;
                uint64_t start = cellKeyCounts[j], end = cellKeyCounts[k];
                if (rangeStart[t] <= end - start && end - start <= rangeEnd[t]) {
                    if (assignKeysBatch(j, start, end, level)) {
                    } else {
                        failedKeyCount += end - start;
                        for (uint64_t i = start; i < end; i++) {
                            failedKeys.push_back(primaryKeyCache[i]);
                            failedValues.push_back(primaryValueCache[i]);
                        }
                    }
                }
            }
            if (failedKeys.empty()) return;
            if (level + 1 == MAXIMUM_LEVEL) return;

            uint64_t failedTotal = failedKeys.size();
            KeyType* nextLevelKeys = new KeyType[failedTotal];
            memset(nextLevelKeys, 0, sizeof(KeyType) * failedTotal);
            ValueType* nextLevelValues = new ValueType[failedTotal];
            memset(nextLevelValues, 0, sizeof(ValueType) * failedTotal);

            for (uint64_t i = 0; i < failedTotal; i++) nextLevelKeys[i] = failedKeys[i];
            for (uint64_t i = 0; i < failedTotal; i++) nextLevelValues[i] = failedValues[i];

            initializeCell(nextLevelKeys, nextLevelValues, failedTotal, level + 1);
            FREE_ARRAY_SPACE(nextLevelKeys);
            FREE_ARRAY_SPACE(nextLevelValues);
        }
    }

    bool assignKeysBatch(uint64_t cellID, uint64_t start, uint64_t end, int level) {
        cellArray[level][cellID] = ((1 << cellSize[level]) - 1);

        uint16_t* cellOffsets = new uint16_t[end - start];

        for (uint64_t i = start; i < end; i++) {
            cellOffsets[i - start] = MurmurHash3_x86_32(&primaryKeyCache[i], sizeof(KeyType), seed2) & cellMask;
            for (uint64_t j = start; j < i; j++)
                if (cellOffsets[j - start] == cellOffsets[i - start]) {
                    FREE_ARRAY_SPACE(cellOffsets);
                    return false;
                }
        }

        for (uint64_t header = 0, maxHeader = (1 << (cellSize[level] - offsetLength[level])); header < maxHeader; header++) {
            bool conflict = false;

            uint64_t availableSlots = (1ull << ((1ull << offsetLength[level]) - 1)) - 1;
            uint32_t baseSlotID = MurmurHash3_x86_32(&cellID, sizeof(uint64_t), seed4 + level + header * 30);
            baseSlotID &= totalSlotMask;

            for (uint64_t i = start; i < end && !conflict; i++) {
                uint64_t slotID = baseSlotID + cellOffsets[i - start];
                if (slotID >= totalSlots) slotID -= totalSlots;
                uint64_t nextSlot1 = (slotID >> 5) + 1;
                if (nextSlot1 == (totalSlots >> 5)) nextSlot1 = 0;
                uint64_t nextSlot2 = nextSlot1 + 1;
                if (nextSlot2 == (totalSlots >> 5)) nextSlot2 = 0;
                uint64_t slotList = ((uint64_t)bitVector[nextSlot2] << 32) | bitVector[nextSlot1];
                slotList = (slotList << (32 - (slotID & 31))) | (bitVector[slotID >> 5] >> (slotID & 31));
                if ((availableSlots & slotList) == 0) {
                    conflict = true;
                }
            }

            if (!conflict) {
                for (uint64_t i = start; i < end; i++) {
                    uint64_t slotID = baseSlotID + cellOffsets[i - start];
                    if (slotID >= totalSlots) slotID -= totalSlots;
                    uint64_t offset = (1ull << (slotID & 31));
                    bitVector[slotID >> 5] |= offset;
                    primaryKeyArray[slotID] = primaryKeyCache[i];
                    primaryValueArray[slotID] = primaryValueCache[i];
                    hierarchyLevels[primaryKeyCache[i]] = level;
                }
                FREE_ARRAY_SPACE(cellOffsets);
                return true;
            }
        }

        FREE_ARRAY_SPACE(cellOffsets);
        return false;
    }

    ValueType queryKeyAtLevel(const KeyType& key, int level, ValueType defaultAnswer = 0) {
        uint64_t cellID = MurmurHash3_x86_32(&key, sizeof(KeyType), seed1) % cellCounts[level];
        
        if (cellArray[level][cellID] == (1 << cellSize[level]) - 1)
            return queryKeyAtLevel(key, level + 1, defaultAnswer);

        uint64_t headerInfo = cellArray[level][cellID] >> offsetLength[level];
        uint32_t baseSlotID = MurmurHash3_x86_32(&cellID, sizeof(uint64_t), seed4 + level + headerInfo * 30);
        baseSlotID &= totalSlotMask;

        uint16_t offsetWithinCell = MurmurHash3_x86_32(&key, sizeof(KeyType), seed2) & cellMask;
        uint64_t slotID = baseSlotID + offsetWithinCell;

        uint64_t realSlotID = slotID + (cellArray[level][cellID] & ((1 << offsetLength[level]) - 1));
        if (realSlotID >= totalSlots) realSlotID -= totalSlots;

        return primaryValueCache[realSlotID];
    }

    uint32_t queryKeyWithInit(const KeyType& key, uint64_t initialCell, uint16_t initialOffset, int level) {
        uint64_t cellID = initialCell;
        uint64_t headerInfo = cellArray[level][cellID] >> offsetLength[level];

        uint32_t baseSlotID = MurmurHash3_x86_32(&cellID, sizeof(uint64_t), seed4 + level + headerInfo * 30);
        baseSlotID &= totalSlotMask;

        uint16_t offsetWithinCell = initialOffset;
        uint64_t slotID = baseSlotID + offsetWithinCell;

        uint64_t realSlotID = slotID + (cellArray[level][cellID] & ((1 << offsetLength[level]) - 1));
        if (realSlotID >= totalSlots) realSlotID -= totalSlots;

        return primaryValueCache[realSlotID];
    }

    bool removeKeyAtLevel(const KeyType& key, int level, ValueType defaultAnswer = 0) {
        uint64_t cellID = MurmurHash3_x86_32(&key, sizeof(KeyType), seed1) % cellCounts[level];
        
        if (cellArray[level][cellID] == (1 << cellSize[level]) - 1)
            return removeKeyAtLevel(key, level + 1, defaultAnswer);

        uint64_t headerInfo = cellArray[level][cellID] >> offsetLength[level];
        uint32_t baseSlotID = MurmurHash3_x86_32(&cellID, sizeof(uint64_t), seed4 + level + headerInfo * 30);
        baseSlotID &= totalSlotMask;

        uint16_t offsetWithinCell = MurmurHash3_x86_32(&key, sizeof(KeyType), seed2) & cellMask;
        uint64_t slotID = baseSlotID + offsetWithinCell;

        uint64_t realSlotID = slotID + (cellArray[level][cellID] & ((1 << offsetLength[level]) - 1));
        if (realSlotID >= totalSlots) realSlotID -= totalSlots;

        primaryValueCache[realSlotID] = ValueType();  
        return true;
    }

    bool removeBitAtLevel(const KeyType& key, int level = 0) {
        uint64_t cellID = MurmurHash3_x86_32(&key, sizeof(KeyType), seed1) % cellCounts[level];
        
        if (cellArray[level][cellID] == (1 << cellSize[level]) - 1)
            return removeBitAtLevel(key, level + 1);

        uint64_t headerInfo = cellArray[level][cellID] >> offsetLength[level];
        uint32_t baseSlotID = MurmurHash3_x86_32(&cellID, sizeof(uint64_t), seed4 + level + headerInfo * 30);
        baseSlotID &= totalSlotMask;

        uint16_t offsetWithinCell = MurmurHash3_x86_32(&key, sizeof(KeyType), seed2) & cellMask;
        uint64_t slotID = baseSlotID + offsetWithinCell;

        uint64_t realSlotID = slotID + (cellArray[level][cellID] & ((1 << offsetLength[level]) - 1));
        if (realSlotID >= totalSlots) realSlotID -= totalSlots;

        bitVector[realSlotID >> 5] &= ~(1ull << (realSlotID & 31));  
        return true;
    }

    uint32_t queryBitAtLevel(const KeyType& key, int level = 0) {
        uint64_t cellIDHashed = MurmurHash3_x86_32(&key, sizeof(KeyType), seed1);
        uint64_t cellID;
        for (int currentLevel = 0; currentLevel < MAXIMUM_LEVEL; ++currentLevel) {
            if (currentLevel > effectiveLevelLimit) {
                return 0;  
            }
            if (currentLevel <= 2) {
                cellID = cellIDHashed % cellCounts[currentLevel];
            }
            if (cellArray[currentLevel][cellID] == (1 << cellSize[currentLevel]) - 1) {
                continue;  
            } else {
                level = currentLevel;
                break;
            }
        }

        uint64_t headerInfo = cellArray[level][cellID] >> offsetLength[level];
        uint32_t baseSlotID = MurmurHash3_x86_32(&cellID, sizeof(uint64_t), seed4 + level + headerInfo * 30);
        baseSlotID &= totalSlotMask;

        uint16_t offsetWithinCell = MurmurHash3_x86_32(&key, sizeof(KeyType), seed2) & cellMask;
        uint64_t slotID = baseSlotID + offsetWithinCell;

        uint64_t realSlotID = slotID + (cellArray[level][cellID] & ((1 << offsetLength[level]) - 1));
        if (realSlotID >= totalSlots) realSlotID -= totalSlots;

        return (bitVector[realSlotID >> 5] >> (realSlotID & 31)) & 1;  
    }

    uint32_t queryBitWithInit(const KeyType& key, uint64_t initialCell, uint16_t initialOffset, int level) {
        uint64_t cellID = initialCell;

        uint64_t headerInfo = cellArray[level][cellID] >> offsetLength[level];
        uint32_t baseSlotID = MurmurHash3_x86_32(&cellID, sizeof(uint64_t), seed4 + level + headerInfo * 30);
        baseSlotID &= totalSlotMask;

        uint16_t offsetWithinCell = initialOffset;
        uint64_t slotID = baseSlotID + offsetWithinCell;

        uint64_t realSlotID = slotID + (cellArray[level][cellID] & ((1 << offsetLength[level]) - 1));
        if (realSlotID >= totalSlots) realSlotID -= totalSlots;

        return (bitVector[realSlotID >> 5] >> (realSlotID & 31)) & 1;  
    }

};


class Segment {
public:
    uint32_t version_lock; 
    ContainerCategory segment_type;
    uint32_t current_key_count;
    uint32_t bitmapSize;
    uint32_t capacity;
    double t_factor;
    double s_factor;
    uint8_t local_depth;
    KeyType* keys;
    ValueType* primaryValueCache;
    BitArray* bitmap;
    // BitArray* availability_filter; 
    // std::map<KeyType, uint8_t> key_to_slot_mapping;
    uint64_t slotCapacity;
    int conflict_counter;
    unsigned char* keyToSlotMapping;
    TransHash* tense_processor;

    Segment(uint32_t initial_capacity = CONTAINER_INIT_SIZE, 
            double _s_factor = DEFAULT_SLT, 
            double _t_factor = DEFAULT_TRANS)
        : capacity(initial_capacity), 
          t_factor(DEFAULT_TRANS), 
          s_factor(_s_factor) {
        Initialize();
    }

    ~Segment() {
        CleanUp();
    }

    void Initialize() {
        keys = new KeyType[capacity];
        primaryValueCache = new ValueType[capacity];
        conflict_counter = 0;
        slotCapacity = capacity << 1; 
        keyToSlotMapping = new unsigned char[slotCapacity];
        std::memset(keyToSlotMapping, 0, slotCapacity * sizeof(unsigned char));
        bitmap = new BitArray();
        bitmap->clearAllBits();
        bitmapSize = bitmap->getSize();
        // availability_filter = new BitArray();
        // availability_filter->clearAllBits();
        segment_type = ContainerCategory::BucketMapContainer;
        local_depth = 1;
        current_key_count = 0;
        version_lock = 0;
    }

    void CleanUp() {
        delete[] keys;
        delete[] primaryValueCache;
        delete bitmap;
        // delete availability_filter;
        delete tense_processor; 
    }

    void RetrieveRandomAccess(uint64_t key, bool& key_found, 
                              uint64_t cell_id = 0, 
                              uint16_t offset = 0, 
                              int level = 0) {
        if (segment_type == ContainerCategory::BucketMapContainer) {
            extractKeyFromBitmap(key, key_found);
        } else {
            uint8_t query_result = tense_processor->queryBitWithInit(key, cell_id, offset, level);
            key_found = (query_result == 1);
        }
    }

    int InsertKey(const KeyType& key, const ValueType& value) {
        if (current_key_count >= capacity * t_factor) {
            return -1; 
        }

        if (!embedKeyInBitmap(key)) {
            return 0; 
        }

        keys[current_key_count] = key;
        primaryValueCache[current_key_count++] = value;


        if (static_cast<double>(current_key_count) >= capacity * t_factor) {
            InitiateTension();
        }

        return 1; 
    }

    void InitiateTension() {
        tense_processor = new TransHash();
        if (segment_type == ContainerCategory::BucketMapContainer) {
            uint64_t total_keys = current_key_count;
            uint64_t total_slots = total_keys / s_factor;
            uint64_t total_cells = 2 * total_keys;
            tense_processor->initialize(keys, primaryValueCache, total_keys, total_slots, total_cells);
            segment_type = ContainerCategory::PerfectHashContainer;
        }
    }

    /* Lock Operations */
    inline void AcquireLock() {
        uint32_t new_version, old_version;
        
        do {
            while (true) {
                old_version = __atomic_load_n(&version_lock, __ATOMIC_ACQUIRE);
                if (!(old_version & lock_pos)) {
                    old_version &= lock_mask;
                    break;
                }
            }
            new_version = old_version | lock_pos;
        } while (!CAS_OPERATION(&version_lock, &old_version, new_version));
    }

    inline bool TryAcquireLock() {
        uint32_t current_version = __atomic_load_n(&version_lock, __ATOMIC_ACQUIRE);
        if (current_version & lock_pos) return false;

        auto old_version = current_version & lock_mask;
        auto new_version = current_version | lock_pos;
        return CAS_OPERATION(&version_lock, &old_version, new_version);
    }

    inline void ReleaseLock() {
        uint32_t current_version = version_lock;
        __atomic_store_n(&version_lock, current_version + 1 - lock_pos, __ATOMIC_RELEASE);
    }

    inline bool IsLockActive(uint32_t& version) {
        version = __atomic_load_n(&version_lock, __ATOMIC_ACQUIRE);
        return (version & lock_pos) != 0;
    }

    inline bool HasLockVersionChanged(uint32_t old_version) {
        auto current_value = __atomic_load_n(&version_lock, __ATOMIC_ACQUIRE);
        return (old_version != current_value);
    }

private:
    uint32_t embedded_missing_count = 0;

    bool embedKeyInBitmap(const KeyType& key, uint8_t startLevel = 0) {
        uint8_t positionShift = 0;

        for (; startLevel < MAXIMUM_LEVEL; ++startLevel) {
            uint32_t initialPosition = (key & seed_of_level[startLevel]) % (bitmapSize << 3);

            if (bitmap->isBitSet(initialPosition)) {
                uint32_t occupancyList = bitmap->extract32BitSlice(initialPosition);
                uint32_t availableSlots = ~occupancyList;

                if (availableSlots == 0) { 
                    if (startLevel == MAXIMUM_LEVEL - 1) return false;
                    continue;  // Move to the next level
                } else {
                    positionShift = __builtin_ctzll(availableSlots);  // Find first available slot
                }
            }

            if (startLevel != 0 || positionShift != 0) {
                uint8_t metaSlot = (startLevel << 5) & positionShift;
                int index = (key & seed_of_level[0]) % slotCapacity;
                keyToSlotMapping[index] = metaSlot;
            }

            bitmap->setBit(initialPosition + positionShift);  // Set the bit to indicate occupancy
            return true;
        }
        return false;
    }

    void extractKeyFromBitmap(const KeyType& key, bool& isKeyFound) {
        uint8_t levelIndicator = 0;
        uint8_t positionShift = 0;

        int index = (key & seed_of_level[0]) % slotCapacity;
        levelIndicator = keyToSlotMapping[index] >> 5;
        positionShift = keyToSlotMapping[index] & 0x1F;

        uint32_t initialPosition = (key & seed_of_level[levelIndicator]) % (bitmapSize << 3);
        isKeyFound = bitmap->isBitSet(initialPosition + positionShift);  // Check if 
    }

};


template<typename KeyType, typename ValueType>
class Container {
public:
    uint32_t version_control_lock;
    uint8_t global_depth;
    uint32_t segment_size;
    Segment*** segment_array;

    Container(uint64_t initial_size = CONTAINER_INIT_SIZE, double s_factor = DEFAULT_SLT)
        : segment_size(initial_size), global_depth(1), version_control_lock(0) {        
        segment_array = new Segment**[1U << global_depth];
        for (uint32_t index = 0; index < (1U << global_depth); ++index) {
            segment_array[index] = new Segment*();
            *(segment_array[index]) = new Segment(segment_size, s_factor);
        }
    }

    void SplitSegment(uint64_t segment_index) {
        uint64_t mask = 1ULL << (global_depth - 1);
        Segment* target_segment = *(segment_array[segment_index]);
        auto zero_segment = new Segment(segment_size);
        auto one_segment = new Segment(segment_size);
        zero_segment->local_depth = target_segment->local_depth + 1;
        one_segment->local_depth = target_segment->local_depth + 1;

        for (uint32_t index = 0; index < target_segment->key_number; index++) {
            KeyType key = target_segment->keys[index];
            ValueType val = target_segment->primaryValueCache[index];
            (key & mask) ? zero_segment->Insert(key, val) : one_segment->Insert(key, val);
        }

        for (uint32_t index = 0; index < (1U << global_depth); ++index) {
            if (*(segment_array[index]) == target_segment) {
                *(segment_array[index]) = (index & mask) ? one_segment : zero_segment;
            } 
        }
    }

    void ExtendContainer() {
        Segment*** old_segment_array = segment_array;
        segment_array = new Segment**[1U << (global_depth + 1)];
        std::memcpy(segment_array, old_segment_array, sizeof(Segment**) * (1U << global_depth));

        for (uint32_t index = 0; index < (1U << global_depth); ++index) {
            segment_array[index + (1U << global_depth)] = new Segment*();
            *(segment_array[index + (1U << global_depth)]) = *(old_segment_array[index]);
        }

        global_depth++;
        delete[] old_segment_array;
    }

    inline void AcquireLock() {
        uint32_t new_version = 0;
        uint32_t old_version = 0;

        do {
            do {
                old_version = __atomic_load_n(&version_control_lock, __ATOMIC_ACQUIRE);
                old_version &= lock_mask;
            } while (old_version & lock_pos);
            new_version = old_version | lock_pos;
        } while (!CAS_OPERATION(&version_control_lock, &old_version, new_version));
    }

    inline bool TryAcquireLock() {
        uint32_t version = __atomic_load_n(&version_control_lock, __ATOMIC_ACQUIRE);
        if (version & lock_pos) return false;

        auto old_version = version & lock_mask;
        auto new_version = version | lock_pos;
        return CAS_OPERATION(&version_control_lock, &old_version, new_version);
    }

    inline void ReleaseLock() {
        __atomic_fetch_add(&version_control_lock, 1 - lock_pos, __ATOMIC_RELEASE);
    }

    inline bool IsLockSet(uint32_t& version) {
        version = __atomic_load_n(&version_control_lock, __ATOMIC_ACQUIRE);
        return (version & lock_pos) != 0;
    }

    inline bool HasLockVersionChanged(uint32_t old_version) {
        return old_version != __atomic_load_n(&version_control_lock, __ATOMIC_ACQUIRE);
    }

    ~Container() {
        for (uint32_t index = 0; index < (1U << global_depth); ++index) {
            delete *(segment_array[index]);
        }
        delete[] segment_array;
    }
};



template<typename KeyType, typename ValueType>
class BitmapContainer {
private:
    std::mutex mtx_container_lock; 
    void InitializeContainerMemory(uint64_t requested_size) {
        key_storage = new KeyType[requested_size];
        value_storage = new ValueType[requested_size];
    }

public:
    uint64_t currently_occupied; 
    uint64_t allocated_capacity; 
    ContainerCategory container_type; 
    KeyType* key_storage; 
    ValueType* value_storage; 

    BitmapContainer(uint64_t initial_size = CONTAINER_INIT_SIZE, const char* filepath = "")
        : currently_occupied(0), allocated_capacity(CONTAINER_INIT_SIZE), container_type(ContainerCategory::BucketMapContainer) {
        if (strcmp(filepath, "") == 0) {
            InitializeContainerMemory(initial_size);
        } else {
            LoadFromFile(filepath);
        }
    }

    bool IsFull() const {
        return allocated_capacity == currently_occupied; 
    }

    void ExpandContainer() {
        Lock<std::mutex> lock(mtx_container_lock);
        ResizeArray<KeyType>(key_storage, allocated_capacity, allocated_capacity << 1);
        ResizeArray<ValueType>(value_storage, allocated_capacity, allocated_capacity << 1);
        allocated_capacity <<= 1;
    }

    void DumpToFile(const char* filepath) const {
        FILE* file_pointer = fopen(filepath, "w");
        char container_type_indicator = (container_type == ContainerCategory::BucketMapContainer) ? 0 : 1;

        fwrite(&container_type_indicator, sizeof(char), 1, file_pointer);
        fwrite(&allocated_capacity, sizeof(uint64_t), 1, file_pointer); 
        fwrite(&currently_occupied, sizeof(uint64_t), 1, file_pointer); 
        fwrite(key_storage, sizeof(KeyType), currently_occupied, file_pointer);
        fwrite(value_storage, sizeof(ValueType), currently_occupied, file_pointer);
        fclose(file_pointer);
    }

    void LoadFromFile(const char* filepath) {
        FILE* file_pointer = fopen(filepath, "r");
        char container_type_indicator;
        fread(&container_type_indicator, sizeof(char), 1, file_pointer);
        
        container_type = (container_type_indicator == 0) ? ContainerCategory::BucketMapContainer : ContainerCategory::PerfectHashContainer;
        fread(&allocated_capacity, sizeof(uint64_t), 1, file_pointer);
        fread(&currently_occupied, sizeof(uint64_t), 1, file_pointer);
        
        InitializeContainerMemory(allocated_capacity);
        fread(key_storage, sizeof(KeyType), currently_occupied, file_pointer);
        fread(value_storage, sizeof(ValueType), currently_occupied, file_pointer);
        fclose(file_pointer);

        allocated_capacity = currently_occupied << 1; 
    }

    void ManualReleaseResources() {
        FreeArrayResources(key_storage);
        FreeArrayResources(value_storage);
    }

    ~BitmapContainer() {
        FreeArrayResources(key_storage);
        FreeArrayResources(value_storage);
    }

    uint64_t Capacity() const {
        return allocated_capacity; 
    }

    void PrintContainerDetails() const {
        std::cout << "Allocated Capacity: " << allocated_capacity << std::endl;
        std::cout << "Currently Occupied: " << currently_occupied << std::endl;
        std::cout << "Container Type: " << (container_type == ContainerCategory::BucketMapContainer ? "BucketMapContainer" : "PerfectHashContainer") << std::endl;
    }
};


template<typename KeyType, typename ValueType>
class BitmapContainerPool {
public:
    BitmapContainerPool(uint64_t pool_size = CONTAINER_POOL_SIZE) 
        : total_capacity(pool_size), available_containers(pool_size) {
        container_pool = new BitmapContainer*[total_capacity];

        for (uint64_t i = 0; i < total_capacity; ++i) {
            container_pool[i] = nullptr;
        }
    }

    BitmapContainer* AcquireContainerByIndex(const uint64_t& index) {
        if (container_pool[index] != nullptr) {
            return container_pool[index];
        } else if (container_pool[index] == nullptr && available_containers != 0) {
            container_pool[index] = new BitmapContainer();
            available_containers--;
            return container_pool[index];
        } else {
            return nullptr;
        }
    }

    BitmapContainer* AcquireSequentialContainer(uint64_t& container_index) {
        if (available_containers == 0) {
            return nullptr; 
        }

        BitmapContainer* candidate_container = container_pool[total_capacity - available_containers];
        if (candidate_container) {
            if (candidate_container->IsFull()) {
                available_containers--;
            } else {
                container_index = total_capacity - available_containers;
                return candidate_container;
            }
        }

        if (available_containers == 0) {
            return nullptr; 
        }

        container_pool[total_capacity - available_containers] = new BitmapContainer();
        container_index = total_capacity - available_containers;
        return container_pool[total_capacity - available_containers];
    }

    void PrintPoolDetails() const {
        for (uint64_t i = 0; i < total_capacity; ++i) {
            if (container_pool[i] != nullptr) {
                std::cout << "Pool[" << i << "]: " << std::endl;
                container_pool[i]->PrintContainerDetails();
            }
        }
    }

    ~BitmapContainerPool() {
        for (uint64_t i = 0; i < total_capacity; ++i) {
            FreeArrayResources(container_pool[i]);
        }
        delete[] container_pool; 
    }

    bool IsEmpty() const {
        return available_containers == 0;
    }

    uint64_t GetPoolSize() const {
        return total_capacity; 
    }

    bool ExpandPool() {
        total_capacity <<= 1; 
        return true; 
    }
    
private:
    BitmapContainer** container_pool; 
    uint64_t total_capacity;
    std::atomic<uint64_t> available_containers; 
    std::mutex mtx_pool_lock; 
};



#endif // !CONTAINER_H