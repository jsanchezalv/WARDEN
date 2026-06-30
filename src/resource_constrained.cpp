#include <Rcpp.h>
#include <queue>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include <cstdint>

using namespace Rcpp;

// Structure to hold patient information in queue with version tracking
struct QueuedPatient {
  int patient_id;
  int64_t priority;
  int insertion_order;
  double queue_start_time;
  uint64_t version;
  int amount_requested;  // units requested when queuing

  QueuedPatient(int id, int64_t prio, int order, double start_time, uint64_t ver, int amt = 1) :
    patient_id(id), priority(prio), insertion_order(order),
    queue_start_time(start_time), version(ver), amount_requested(amt) {}
};

// Comparator: higher priority first, then FIFO or LIFO within same priority
struct QueueComparator {
  bool is_lifo;
  QueueComparator(bool lifo = false) : is_lifo(lifo) {}
  bool operator()(const QueuedPatient& a, const QueuedPatient& b) const {
    if (a.priority != b.priority) return a.priority < b.priority;
    if (is_lifo) return a.insertion_order < b.insertion_order;  // LIFO: higher order wins
    return a.insertion_order > b.insertion_order;               // FIFO: lower order wins
  }
};

// Structure to hold patient information when using resource
struct UsingPatient {
  int patient_id;
  double start_time;
  int amount;  // units occupied by this patient

  UsingPatient(int id, double time, int amt = 1) : patient_id(id), start_time(time), amount(amt) {}
};

class DiscreteResource {
private:
  int total_capacity;
  int current_total_used;  // maintained incrementally for O(1) n_free()
  bool is_lifo;
  int max_queue_capacity;  // -1 = unlimited
  bool allow_multiple_queue;
  int next_insertion_order;
  int64_t current_max_priority;
  int total_entries_ever_added;
  int current_valid_entries;
  int operations_since_cleanup;
  uint64_t next_version;

  std::vector<UsingPatient> patients_using;
  std::priority_queue<QueuedPatient, std::vector<QueuedPatient>, QueueComparator> patient_queue;

  std::unordered_map<int, int> patient_queue_count;
  std::unordered_map<int, std::vector<double>> patient_queue_start_times;
  std::unordered_map<int, uint64_t> patient_current_version;

  // Statistics
  std::unordered_set<int> unique_patients_blocked;
  std::unordered_set<int> unique_patients_queued;
  std::unordered_map<int, double> patient_last_queue_wait;  // patient_id → last wait

  static const int MAX_QUEUE_ENTRIES_PER_PATIENT = 1000;
  static const int CLEANUP_FREQUENCY = 10000;

  bool is_entry_valid(const QueuedPatient& entry) const {
    auto it = patient_current_version.find(entry.patient_id);
    return (it != patient_current_version.end()) && (entry.version == it->second);
  }

  void cleanup_queue_top() {
    while (!patient_queue.empty() && !is_entry_valid(patient_queue.top())) {
      patient_queue.pop();
    }
  }

  void full_queue_cleanup() {
    std::vector<QueuedPatient> valid_patients;

    while (!patient_queue.empty()) {
      QueuedPatient patient = patient_queue.top();
      patient_queue.pop();
      if (is_entry_valid(patient)) valid_patients.push_back(patient);
    }

    for (const auto& patient : valid_patients) patient_queue.push(patient);

    total_entries_ever_added = static_cast<int>(valid_patients.size());
    current_valid_entries = static_cast<int>(valid_patients.size());
  }

  void check_and_cleanup() {
    operations_since_cleanup++;

    bool periodic_cleanup = (operations_since_cleanup >= CLEANUP_FREQUENCY);
    bool threshold_cleanup = (total_entries_ever_added - current_valid_entries >
                                static_cast<int>(current_valid_entries * 2));

    if (periodic_cleanup || threshold_cleanup) {
      full_queue_cleanup();
      operations_since_cleanup = 0;
    }
  }

public:
  DiscreteResource(int n, bool lifo = false, int max_queue = -1,
                   bool allow_multi_queue = true) :
    total_capacity(n),
    current_total_used(0),
    is_lifo(lifo),
    max_queue_capacity(max_queue),
    allow_multiple_queue(allow_multi_queue),
    next_insertion_order(0),
    current_max_priority(1),
    total_entries_ever_added(0),
    current_valid_entries(0),
    operations_since_cleanup(0),
    next_version(1),
    patient_queue(QueueComparator(lifo)) {
    if (n < 0) stop("Resource capacity must be >= 0");
    if (n > 0) patients_using.reserve(n);
  }

  int size() const { return total_capacity; }
  int queue_size() const { return current_valid_entries; }
  int n_free() const { return total_capacity - current_total_used; }

  std::vector<int> patients_using_ids() const {
    std::vector<int> ids;
    ids.reserve(patients_using.size());
    for (const auto& p : patients_using) ids.push_back(p.patient_id);
    return ids;
  }

  std::vector<double> patients_using_times() const {
    std::vector<double> times;
    times.reserve(patients_using.size());
    for (const auto& p : patients_using) times.push_back(p.start_time);
    return times;
  }

  bool is_patient_in_queue(int patient_id) const {
    return patient_current_version.find(patient_id) != patient_current_version.end();
  }

  bool is_patient_using(int patient_id) const {
    return std::find_if(patients_using.begin(), patients_using.end(),
                        [patient_id](const UsingPatient& p) {
                          return p.patient_id == patient_id;
                        }) != patients_using.end();
  }

  double get_patient_using_start_time(int patient_id) const {
    auto it = std::find_if(patients_using.begin(), patients_using.end(),
                           [patient_id](const UsingPatient& p) {
                             return p.patient_id == patient_id;
                           });
    if (it == patients_using.end()) return NA_REAL;
    return it->start_time;
  }

  bool can_patient_acquire(int patient_id, int amount) {
    if (n_free() < amount) return false;
    if (current_valid_entries == 0) return true;
    cleanup_queue_top();
    if (patient_queue.empty()) return true;
    const QueuedPatient& top = patient_queue.top();
    return (top.patient_id == patient_id && is_entry_valid(top));
  }

  bool can_patient_queue(int patient_id) {
    if (max_queue_capacity >= 0 && current_valid_entries >= max_queue_capacity) return false;
    if (!allow_multiple_queue &&
        patient_queue_count.find(patient_id) != patient_queue_count.end() &&
        patient_queue_count.at(patient_id) > 0) return false;
    return true;
  }

  // Returns 1 (acquired), 0 (queued), -1 (rejected)
  int attempt_block(int patient_id, int64_t priority, double start_time, int amount = 1) {
    if (patient_queue_count.find(patient_id) != patient_queue_count.end() &&
        patient_queue_count[patient_id] >= MAX_QUEUE_ENTRIES_PER_PATIENT) {
      stop("Patient exceeds maximum queue entries limit");
    }

    // Direct acquire: no queue and enough capacity
    if (current_valid_entries == 0 && n_free() >= amount) {
      patients_using.emplace_back(patient_id, start_time, amount);
      current_total_used += amount;
      unique_patients_blocked.insert(patient_id);
      return 1;
    }

    // Dequeue-and-acquire: patient is first in line and capacity covers their queued amount.
    // Capacity is checked against the QUEUED amount (immutable), not the retry call's amount.
    if (current_valid_entries > 0) {
      cleanup_queue_top();

      if (!patient_queue.empty()) {
        QueuedPatient next_in_line = patient_queue.top();
        if (next_in_line.patient_id == patient_id && is_entry_valid(next_in_line) &&
            n_free() >= next_in_line.amount_requested) {

          // Compute and store wait time before erasing queue data
          double queue_entry_time = 0.0;
          if (patient_queue_start_times.find(patient_id) != patient_queue_start_times.end() &&
              !patient_queue_start_times[patient_id].empty()) {
            queue_entry_time = patient_queue_start_times[patient_id].front();
          }
          patient_last_queue_wait[patient_id] = start_time - queue_entry_time;

          patient_queue.pop();
          patients_using.emplace_back(patient_id, start_time, next_in_line.amount_requested);
          current_total_used += next_in_line.amount_requested;

          patient_queue_count[patient_id]--;
          if (patient_queue_count[patient_id] == 0) {
            patient_queue_count.erase(patient_id);
            patient_queue_start_times.erase(patient_id);
            patient_current_version.erase(patient_id);
          } else {
            patient_queue_start_times[patient_id].erase(patient_queue_start_times[patient_id].begin());
          }
          current_valid_entries--;
          unique_patients_blocked.insert(patient_id);
          return 1;
        }
      }
    }

    // Queue or reject
    // Reject if allow_multiple_queue = FALSE and patient already has a queue entry
    if (!allow_multiple_queue &&
        patient_queue_count.find(patient_id) != patient_queue_count.end() &&
        patient_queue_count[patient_id] > 0) {
      return -1;
    }

    if (max_queue_capacity >= 0 && current_valid_entries >= max_queue_capacity) {
      return -1;  // queue full — reject
    }

    uint64_t patient_version;
    if (patient_current_version.find(patient_id) == patient_current_version.end()) {
      patient_version = next_version++;
      patient_current_version[patient_id] = patient_version;
      patient_queue_count[patient_id] = 0;
      patient_queue_start_times[patient_id].reserve(10);
    } else {
      patient_version = patient_current_version[patient_id];
    }

    patient_queue.emplace(patient_id, priority, next_insertion_order++, start_time, patient_version, amount);

    patient_queue_count[patient_id]++;
    patient_queue_start_times[patient_id].push_back(start_time);

    current_max_priority = std::max(current_max_priority, priority);

    total_entries_ever_added++;
    current_valid_entries++;

    unique_patients_queued.insert(patient_id);

    check_and_cleanup();
    return 0;
  }

  void attempt_free(int patient_id, bool remove_all = false, int amount = NA_INTEGER) {
    // NA_INTEGER sentinel means amount = 1 (matches release(amount = NULL) semantics)
    if (amount == NA_INTEGER) amount = 1;

    bool found_in_using = false;

    if (remove_all) {
      int freed = 0;
      for (const auto& p : patients_using) {
        if (p.patient_id == patient_id) freed += p.amount;
      }
      auto original_size = patients_using.size();
      patients_using.erase(
        std::remove_if(patients_using.begin(), patients_using.end(),
                       [patient_id](const UsingPatient& p) { return p.patient_id == patient_id; }),
        patients_using.end());
      found_in_using = original_size != patients_using.size();
      if (found_in_using) current_total_used -= freed;
    } else {
      // Indivisible model: find OLDEST entry with matching amount
      auto it = std::find_if(patients_using.begin(), patients_using.end(),
                             [patient_id, amount](const UsingPatient& p) {
                               return p.patient_id == patient_id && p.amount == amount;
                             });
      if (it != patients_using.end()) {
        current_total_used -= it->amount;
        patients_using.erase(it);
        found_in_using = true;
      } else {
        // Check if patient is using with a different amount → indivisible mismatch error
        bool patient_is_using = std::any_of(patients_using.begin(), patients_using.end(),
                                            [patient_id](const UsingPatient& p) {
                                              return p.patient_id == patient_id;
                                            });
        if (patient_is_using) {
          stop("release() amount mismatch: patient is using the resource but not with the specified amount (indivisible units)");
        }
      }
    }

    if (found_in_using) return;

    // Patient not in using list: target queue entries
    if (patient_queue_count.find(patient_id) != patient_queue_count.end() &&
        patient_queue_count[patient_id] > 0) {

      if (remove_all) {
        current_valid_entries -= patient_queue_count[patient_id];
        patient_queue_count.erase(patient_id);
        patient_queue_start_times.erase(patient_id);
        patient_current_version.erase(patient_id);
      } else {
        patient_queue_count[patient_id]--;
        current_valid_entries--;

        if (!patient_queue_start_times[patient_id].empty()) {
          patient_queue_start_times[patient_id].erase(patient_queue_start_times[patient_id].begin());
        }

        if (patient_queue_count[patient_id] == 0) {
          patient_queue_count.erase(patient_id);
          patient_queue_start_times.erase(patient_id);
          patient_current_version.erase(patient_id);
        }
      }

      check_and_cleanup();
    }
  }

  void purge_from_queue(int patient_id) {
    auto it = patient_queue_count.find(patient_id);
    if (it != patient_queue_count.end() && it->second > 0) {
      current_valid_entries -= it->second;
      patient_queue_count.erase(it);
      patient_queue_start_times.erase(patient_id);
      patient_current_version.erase(patient_id);
      check_and_cleanup();
    }
  }

  bool release_full(int patient_id, int amount = NA_INTEGER, bool purge_queue = true) {
    bool was_using = false;

    if (amount == NA_INTEGER) {
      // Release ALL usage entries for this patient
      int freed = 0;
      for (const auto& p : patients_using) {
        if (p.patient_id == patient_id) freed += p.amount;
      }
      auto original_size = patients_using.size();
      patients_using.erase(
        std::remove_if(patients_using.begin(), patients_using.end(),
                       [patient_id](const UsingPatient& p) { return p.patient_id == patient_id; }),
        patients_using.end());
      was_using = original_size != patients_using.size();
      if (was_using) current_total_used -= freed;
    } else {
      // Release ONE matching entry (indivisible)
      auto it = std::find_if(patients_using.begin(), patients_using.end(),
                             [patient_id, amount](const UsingPatient& p) {
                               return p.patient_id == patient_id && p.amount == amount;
                             });
      if (it != patients_using.end()) {
        current_total_used -= it->amount;
        patients_using.erase(it);
        was_using = true;
      } else {
        bool patient_is_using = std::any_of(patients_using.begin(), patients_using.end(),
                                            [patient_id](const UsingPatient& p) {
                                              return p.patient_id == patient_id;
                                            });
        if (patient_is_using) {
          stop("release_all() amount mismatch: patient is using the resource but not with the specified amount (indivisible units)");
        }
      }
    }

    if (purge_queue) {
      purge_from_queue(patient_id);
    }

    return was_using;
  }

  void attempt_free_if_using(int patient_id, bool remove_all = false, int amount = NA_INTEGER) {
    // NA_INTEGER → release ALL usage entries for the patient
    if (amount == NA_INTEGER) {
      int freed = 0;
      for (const auto& p : patients_using) {
        if (p.patient_id == patient_id) freed += p.amount;
      }
      auto original_size = patients_using.size();
      patients_using.erase(
        std::remove_if(patients_using.begin(), patients_using.end(),
                       [patient_id](const UsingPatient& p) { return p.patient_id == patient_id; }),
        patients_using.end());
      if (original_size != patients_using.size()) current_total_used -= freed;
      return;
    }

    if (remove_all) {
      int freed = 0;
      for (const auto& p : patients_using) {
        if (p.patient_id == patient_id) freed += p.amount;
      }
      auto original_size = patients_using.size();
      patients_using.erase(
        std::remove_if(patients_using.begin(), patients_using.end(),
                       [patient_id](const UsingPatient& p) { return p.patient_id == patient_id; }),
        patients_using.end());
      if (original_size != patients_using.size()) current_total_used -= freed;
    } else {
      // Indivisible: find oldest entry matching exact amount
      auto it = std::find_if(patients_using.begin(), patients_using.end(),
                             [patient_id, amount](const UsingPatient& p) {
                               return p.patient_id == patient_id && p.amount == amount;
                             });
      if (it != patients_using.end()) {
        current_total_used -= it->amount;
        patients_using.erase(it);
      } else {
        bool patient_is_using = std::any_of(patients_using.begin(), patients_using.end(),
                                            [patient_id](const UsingPatient& p) {
                                              return p.patient_id == patient_id;
                                            });
        if (patient_is_using) {
          stop("release_all_if_using() amount mismatch: patient is using the resource but not with the specified amount (indivisible units)");
        }
      }
    }
  }

  std::vector<int> next_patient_in_line(int n = 1) {
    std::vector<int> result;
    QueueComparator cmp(is_lifo);
    std::priority_queue<QueuedPatient, std::vector<QueuedPatient>, QueueComparator> temp_queue(cmp);

    int items_to_check = std::min(n * 10, static_cast<int>(patient_queue.size()));
    if (items_to_check < 50) items_to_check = std::min(50, static_cast<int>(patient_queue.size()));

    std::priority_queue<QueuedPatient, std::vector<QueuedPatient>, QueueComparator> original_queue = patient_queue;

    for (int i = 0; i < items_to_check && !original_queue.empty(); ++i) {
      temp_queue.push(original_queue.top());
      original_queue.pop();
    }

    int count = std::min(n, current_valid_entries);
    int checked = 0;
    const int max_check = static_cast<int>(patient_queue.size());

    while (static_cast<int>(result.size()) < count && !temp_queue.empty() && checked < max_check) {
      QueuedPatient patient = temp_queue.top();
      temp_queue.pop();
      checked++;

      if (is_entry_valid(patient)) result.push_back(patient.patient_id);

      if (temp_queue.empty() && static_cast<int>(result.size()) < count && !original_queue.empty()) {
        int additional_items = std::min(50, static_cast<int>(original_queue.size()));
        for (int i = 0; i < additional_items && !original_queue.empty(); ++i) {
          temp_queue.push(original_queue.top());
          original_queue.pop();
        }
      }
    }

    return result;
  }

  std::vector<int> queue_priorities() {
    std::vector<int> result;
    QueueComparator cmp2(is_lifo);
    std::priority_queue<QueuedPatient, std::vector<QueuedPatient>, QueueComparator> temp_queue(cmp2);

    int items_to_check = std::min(current_valid_entries * 10, static_cast<int>(patient_queue.size()));
    if (items_to_check < 50) items_to_check = std::min(50, static_cast<int>(patient_queue.size()));

    std::priority_queue<QueuedPatient, std::vector<QueuedPatient>, QueueComparator> original_queue = patient_queue;

    for (int i = 0; i < items_to_check && !original_queue.empty(); ++i) {
      temp_queue.push(original_queue.top());
      original_queue.pop();
    }

    int count = current_valid_entries;
    int checked = 0;
    const int max_check = static_cast<int>(patient_queue.size());

    while (static_cast<int>(result.size()) < count && !temp_queue.empty() && checked < max_check) {
      QueuedPatient patient = temp_queue.top();
      temp_queue.pop();
      checked++;

      if (is_entry_valid(patient)) result.push_back(static_cast<int>(patient.priority));

      if (temp_queue.empty() && static_cast<int>(result.size()) < count && !original_queue.empty()) {
        int additional_items = std::min(50, static_cast<int>(original_queue.size()));
        for (int i = 0; i < additional_items && !original_queue.empty(); ++i) {
          temp_queue.push(original_queue.top());
          original_queue.pop();
        }
      }
    }

    return result;
  }

  std::vector<double> queue_start_times() {
    std::vector<double> result;
    QueueComparator cmp3(is_lifo);
    std::priority_queue<QueuedPatient, std::vector<QueuedPatient>, QueueComparator> temp_queue(cmp3);
    std::priority_queue<QueuedPatient, std::vector<QueuedPatient>, QueueComparator> original_queue = patient_queue;

    int items_to_check = std::min(current_valid_entries * 2, static_cast<int>(patient_queue.size()));
    for (int i = 0; i < items_to_check && !original_queue.empty(); ++i) {
      temp_queue.push(original_queue.top());
      original_queue.pop();
    }

    std::unordered_map<int, int> patient_index;

    while (static_cast<int>(result.size()) < current_valid_entries && !temp_queue.empty()) {
      QueuedPatient patient = temp_queue.top();
      temp_queue.pop();

      if (is_entry_valid(patient)) {
        int& index = patient_index[patient.patient_id];
        if (patient_queue_start_times.find(patient.patient_id) != patient_queue_start_times.end() &&
            index < static_cast<int>(patient_queue_start_times.at(patient.patient_id).size())) {
          result.push_back(patient_queue_start_times.at(patient.patient_id)[index]);
        } else {
          result.push_back(patient.queue_start_time);
        }
        index++;
      }
    }

    return result;
  }

  void modify_priority(int patient_id, int64_t new_priority) {
    if (patient_queue_count.find(patient_id) == patient_queue_count.end() ||
        patient_queue_count[patient_id] == 0) {
      return;
    }

    int count = patient_queue_count[patient_id];
    std::vector<double> start_times = patient_queue_start_times[patient_id];

    // Retrieve stored amounts before invalidating
    std::vector<int> amounts;
    {
      std::priority_queue<QueuedPatient, std::vector<QueuedPatient>, QueueComparator> tmp = patient_queue;
      while (!tmp.empty()) {
        const QueuedPatient& top = tmp.top();
        if (top.patient_id == patient_id && is_entry_valid(top)) {
          amounts.push_back(top.amount_requested);
        }
        tmp.pop();
      }
    }
    while (static_cast<int>(amounts.size()) < count) amounts.push_back(1);

    patient_current_version[patient_id] = next_version++;
    uint64_t new_version = patient_current_version[patient_id];

    for (int i = 0; i < count; ++i) {
      patient_queue.emplace(patient_id, new_priority, next_insertion_order++,
                            start_times[i], new_version, amounts[i]);
      total_entries_ever_added++;
    }

    check_and_cleanup();
  }

  void move_to_front(int patient_id) {
    if (patient_queue_count.find(patient_id) == patient_queue_count.end() ||
        patient_queue_count[patient_id] == 0) {
      return;
    }

    int count = patient_queue_count[patient_id];
    std::vector<double> start_times = patient_queue_start_times[patient_id];

    // Collect amounts from existing entries
    std::vector<int> amounts;
    {
      std::priority_queue<QueuedPatient, std::vector<QueuedPatient>, QueueComparator> tmp = patient_queue;
      while (!tmp.empty()) {
        const QueuedPatient& top = tmp.top();
        if (top.patient_id == patient_id && is_entry_valid(top)) {
          amounts.push_back(top.amount_requested);
        }
        tmp.pop();
      }
    }
    while (static_cast<int>(amounts.size()) < count) amounts.push_back(1);

    // Use a priority above current max to jump to front
    int64_t elevated_priority = current_max_priority + 1;
    current_max_priority = elevated_priority;

    // Bump version to invalidate old entries
    patient_current_version[patient_id] = next_version++;
    uint64_t new_version = patient_current_version[patient_id];

    for (int j = 0; j < count; ++j) {
      patient_queue.emplace(patient_id, elevated_priority, next_insertion_order++,
                            start_times[j], new_version, amounts[j]);
      total_entries_ever_added++;
    }

    check_and_cleanup();
  }

  void add_resource(int n_to_add) {
    if (n_to_add <= 0) stop("n_to_add must be positive");
    total_capacity += n_to_add;
  }

  void remove_resource(int n_to_remove, double current_time) {
    if (n_to_remove <= 0) stop("n_to_remove must be positive");
    if (n_to_remove > total_capacity) stop("Cannot remove more resources than available");

    total_capacity -= n_to_remove;

    while (current_total_used > total_capacity) {
      UsingPatient patient = patients_using.back();
      patients_using.pop_back();
      current_total_used -= patient.amount;

      int64_t new_priority = current_max_priority + 1;
      current_max_priority = new_priority;

      patient_current_version[patient.patient_id] = next_version++;
      uint64_t patient_version = patient_current_version[patient.patient_id];

      if (patient_queue_count.find(patient.patient_id) == patient_queue_count.end()) {
        patient_queue_count[patient.patient_id] = 0;
        patient_queue_start_times[patient.patient_id].reserve(10);
      }

      patient_queue.emplace(patient.patient_id, new_priority, next_insertion_order++,
                            current_time, patient_version, patient.amount);

      patient_queue_count[patient.patient_id]++;
      patient_queue_start_times[patient.patient_id].push_back(current_time);

      total_entries_ever_added++;
      current_valid_entries++;
    }

    check_and_cleanup();
  }

  // Settings accessors (used by clone)
  bool get_is_lifo() const { return is_lifo; }
  int get_max_queue_capacity() const { return max_queue_capacity; }
  bool get_allow_multiple_queue() const { return allow_multiple_queue; }

  // Statistics accessors
  double queue_wait_time(int patient_id) const {
    auto it = patient_last_queue_wait.find(patient_id);
    if (it == patient_last_queue_wait.end()) return NA_REAL;
    return it->second;
  }

  // Elapsed wait: current_time - entry if still in queue, stored final wait if dequeued, NA if never queued.
  double queue_elapsed_time(int patient_id, double current_time) const {
    auto it_queue = patient_queue_start_times.find(patient_id);
    if (it_queue != patient_queue_start_times.end() && !it_queue->second.empty()) {
      return current_time - it_queue->second.front();
    }
    auto it_wait = patient_last_queue_wait.find(patient_id);
    if (it_wait != patient_last_queue_wait.end()) return it_wait->second;
    return NA_REAL;
  }

  int had_to_queue(int patient_id) const {
    return unique_patients_queued.count(patient_id) > 0 ? 1 : 0;
  }

  int total_patients_blocked() const {
    return static_cast<int>(unique_patients_blocked.size());
  }

  int total_patients_queued() const {
    return static_cast<int>(unique_patients_queued.size());
  }
};

// ─────────────────────────────────────────────────────────────────────────────
// Rcpp exports
// ─────────────────────────────────────────────────────────────────────────────

// [[Rcpp::export]]
SEXP create_discrete_resource_cpp(int n, bool lifo = false, int max_queue_capacity = -1,
                                   bool allow_multiple_queue = true) {
  DiscreteResource* ptr = new DiscreteResource(n, lifo, max_queue_capacity, allow_multiple_queue);
  XPtr<DiscreteResource> xptr(ptr, true);
  return xptr;
}

void validate_xptr(SEXP xptr) {
  if (TYPEOF(xptr) != EXTPTRSXP) stop("Invalid external pointer");
}

// [[Rcpp::export]]
int discrete_resource_size_cpp(SEXP xptr) {
  validate_xptr(xptr);
  return XPtr<DiscreteResource>(xptr)->size();
}

// [[Rcpp::export]]
int discrete_resource_queue_size_cpp(SEXP xptr) {
  validate_xptr(xptr);
  return XPtr<DiscreteResource>(xptr)->queue_size();
}

// [[Rcpp::export]]
int discrete_resource_n_free_cpp(SEXP xptr) {
  validate_xptr(xptr);
  return XPtr<DiscreteResource>(xptr)->n_free();
}

// [[Rcpp::export]]
IntegerVector discrete_resource_patients_using_cpp(SEXP xptr) {
  validate_xptr(xptr);
  return wrap(XPtr<DiscreteResource>(xptr)->patients_using_ids());
}

// [[Rcpp::export]]
NumericVector discrete_resource_patients_using_times_cpp(SEXP xptr) {
  validate_xptr(xptr);
  return wrap(XPtr<DiscreteResource>(xptr)->patients_using_times());
}

// [[Rcpp::export]]
bool discrete_resource_is_patient_in_queue_cpp(SEXP xptr, int patient_id) {
  validate_xptr(xptr);
  return XPtr<DiscreteResource>(xptr)->is_patient_in_queue(patient_id);
}

// [[Rcpp::export]]
bool discrete_resource_is_patient_using_cpp(SEXP xptr, int patient_id) {
  validate_xptr(xptr);
  return XPtr<DiscreteResource>(xptr)->is_patient_using(patient_id);
}

// Returns 1 (acquired), 0 (queued), -1 (rejected)
// [[Rcpp::export]]
int discrete_resource_attempt_block_cpp(SEXP xptr, int patient_id, int priority,
                                         double start_time, int amount = 1) {
  validate_xptr(xptr);
  return XPtr<DiscreteResource>(xptr)->attempt_block(patient_id, priority, start_time, amount);
}

// [[Rcpp::export]]
void discrete_resource_attempt_free_cpp(SEXP xptr, int patient_id,
                                         bool remove_all = false, int amount = NA_INTEGER) {
  validate_xptr(xptr);
  XPtr<DiscreteResource>(xptr)->attempt_free(patient_id, remove_all, amount);
}

// [[Rcpp::export]]
void discrete_resource_attempt_free_if_using_cpp(SEXP xptr, int patient_id,
                                                   bool remove_all = false,
                                                   int amount = NA_INTEGER) {
  validate_xptr(xptr);
  XPtr<DiscreteResource>(xptr)->attempt_free_if_using(patient_id, remove_all, amount);
}

// [[Rcpp::export]]
IntegerVector discrete_resource_next_patient_in_line_cpp(SEXP xptr, int n = 1) {
  validate_xptr(xptr);
  return wrap(XPtr<DiscreteResource>(xptr)->next_patient_in_line(n));
}

// [[Rcpp::export]]
IntegerVector discrete_resource_queue_priorities_cpp(SEXP xptr) {
  validate_xptr(xptr);
  return wrap(XPtr<DiscreteResource>(xptr)->queue_priorities());
}

// [[Rcpp::export]]
NumericVector discrete_resource_queue_start_times_cpp(SEXP xptr) {
  validate_xptr(xptr);
  return wrap(XPtr<DiscreteResource>(xptr)->queue_start_times());
}

// [[Rcpp::export]]
void discrete_resource_modify_priority_cpp(SEXP xptr, int patient_id, int new_priority) {
  validate_xptr(xptr);
  XPtr<DiscreteResource>(xptr)->modify_priority(patient_id, new_priority);
}

// [[Rcpp::export]]
void discrete_resource_add_resource_cpp(SEXP xptr, int n_to_add) {
  validate_xptr(xptr);
  XPtr<DiscreteResource>(xptr)->add_resource(n_to_add);
}

// [[Rcpp::export]]
void discrete_resource_remove_resource_cpp(SEXP xptr, int n_to_remove, double current_time) {
  validate_xptr(xptr);
  XPtr<DiscreteResource>(xptr)->remove_resource(n_to_remove, current_time);
}

// ── New statistics accessors ─────────────────────────────────────────────────

// [[Rcpp::export]]
double discrete_resource_queue_wait_time_cpp(SEXP xptr, int patient_id) {
  validate_xptr(xptr);
  return XPtr<DiscreteResource>(xptr)->queue_wait_time(patient_id);
}

// [[Rcpp::export]]
double discrete_resource_queue_elapsed_time_cpp(SEXP xptr, int patient_id, double current_time) {
  validate_xptr(xptr);
  return XPtr<DiscreteResource>(xptr)->queue_elapsed_time(patient_id, current_time);
}

// [[Rcpp::export]]
int discrete_resource_had_to_queue_cpp(SEXP xptr, int patient_id) {
  validate_xptr(xptr);
  return XPtr<DiscreteResource>(xptr)->had_to_queue(patient_id);
}

// [[Rcpp::export]]
double discrete_resource_get_patient_using_start_time_cpp(SEXP xptr, int patient_id) {
  validate_xptr(xptr);
  return XPtr<DiscreteResource>(xptr)->get_patient_using_start_time(patient_id);
}

// [[Rcpp::export]]
int discrete_resource_total_patients_blocked_cpp(SEXP xptr) {
  validate_xptr(xptr);
  return XPtr<DiscreteResource>(xptr)->total_patients_blocked();
}

// [[Rcpp::export]]
int discrete_resource_total_patients_queued_cpp(SEXP xptr) {
  validate_xptr(xptr);
  return XPtr<DiscreteResource>(xptr)->total_patients_queued();
}

// ── Batch seize ──────────────────────────────────────────────────────────────

// [[Rcpp::export]]
IntegerVector discrete_resource_batch_seize_cpp(SEXP xptr, IntegerVector patient_ids,
                                                  int priority, double start_time,
                                                  int amount_each = 1) {
  validate_xptr(xptr);
  XPtr<DiscreteResource> ptr(xptr);
  IntegerVector results(patient_ids.size());
  for (int i = 0; i < patient_ids.size(); ++i) {
    results[i] = ptr->attempt_block(patient_ids[i], priority, start_time, amount_each);
  }
  return results;
}

// ── Multi-resource seize / release ───────────────────────────────────────────

// [[Rcpp::export]]
int discrete_resource_seize_all_cpp(List resource_xptrs, int patient_id,
                                     IntegerVector priorities, double start_time,
                                     IntegerVector amounts, int policy,
                                     bool force_unblock = false,
                                     bool accum_queue = true) {
  int n = resource_xptrs.size();

  if (policy == 0) {  // all_or_none
    // Phase 1: check what each resource would do for this patient
    std::vector<bool> can_acquire(n);
    for (int i = 0; i < n; ++i) {
      XPtr<DiscreteResource> ptr(VECTOR_ELT(resource_xptrs, i));
      can_acquire[i] = ptr->can_patient_acquire(patient_id, amounts[i]);
    }

    // All can acquire: acquire all atomically
    bool all_ok = true;
    for (int i = 0; i < n; ++i) if (!can_acquire[i]) { all_ok = false; break; }
    if (all_ok) {
      for (int i = 0; i < n; ++i) {
        XPtr<DiscreteResource> ptr(VECTOR_ELT(resource_xptrs, i));
        ptr->attempt_block(patient_id, priorities[i], start_time, amounts[i]);
      }
      return 1;
    }

    // Some cannot acquire: check for deadlock when force_unblock = TRUE
    // Deadlock: all resources have sufficient capacity but patient is not first on some queues
    if (force_unblock) {
      bool all_have_capacity = true;
      for (int i = 0; i < n; ++i) {
        XPtr<DiscreteResource> ptr(VECTOR_ELT(resource_xptrs, i));
        if (ptr->n_free() < amounts[i]) { all_have_capacity = false; break; }
      }
      if (all_have_capacity) {
        // Verify it is a true deadlock: patient must already be queued on every
        // resource where they cannot directly acquire. If they are absent from any
        // such queue, this is a first-call scenario (not a deadlock), so fall through.
        bool is_true_deadlock = true;
        for (int i = 0; i < n; ++i) {
          if (!can_acquire[i]) {
            XPtr<DiscreteResource> ptr(VECTOR_ELT(resource_xptrs, i));
            if (!ptr->is_patient_in_queue(patient_id)) { is_true_deadlock = false; break; }
          }
        }
        if (is_true_deadlock) {
          // Move patient to front on all resources where they are not already first
          for (int i = 0; i < n; ++i) {
            if (!can_acquire[i]) {
              XPtr<DiscreteResource> ptr(VECTOR_ELT(resource_xptrs, i));
              ptr->move_to_front(patient_id);
            }
          }
          // Now acquire all
          for (int i = 0; i < n; ++i) {
            XPtr<DiscreteResource> ptr(VECTOR_ELT(resource_xptrs, i));
            ptr->attempt_block(patient_id, priorities[i], start_time, amounts[i]);
          }
          return 1;
        }
        // Not a true deadlock: patient not yet queued on some resources, fall through
      }
      // Capacity insufficient even with force_unblock: fall through to normal queuing
    }

    // Queue on ALL bottleneck resources
    // Determine which resources need a new queue entry
    std::vector<bool> needs_new_entry(n, false);
    for (int i = 0; i < n; ++i) {
      if (!can_acquire[i]) {
        XPtr<DiscreteResource> ptr(VECTOR_ELT(resource_xptrs, i));
        bool already_queued = ptr->is_patient_in_queue(patient_id);
        if (already_queued) {
          if (!accum_queue || !ptr->get_allow_multiple_queue()) {
            needs_new_entry[i] = false;
          } else {
            needs_new_entry[i] = ptr->can_patient_queue(patient_id);
          }
        } else {
          if (!ptr->can_patient_queue(patient_id)) return -1;
          needs_new_entry[i] = true;
        }
      }
    }

    for (int i = 0; i < n; ++i) {
      if (needs_new_entry[i]) {
        XPtr<DiscreteResource> ptr(VECTOR_ELT(resource_xptrs, i));
        ptr->attempt_block(patient_id, priorities[i], start_time, amounts[i]);
      }
    }
    return 0;

  } else {  // sequential (unchanged)
    for (int i = 0; i < n; ++i) {
      XPtr<DiscreteResource> ptr(VECTOR_ELT(resource_xptrs, i));
      int result = ptr->attempt_block(patient_id, priorities[i], start_time, amounts[i]);
      if (result != 1) return 0;
    }
    return 1;
  }
}

// [[Rcpp::export]]
LogicalVector discrete_resource_release_all_cpp(List resource_xptrs, int patient_id,
                                                 IntegerVector amounts,
                                                 bool purge_queues = true) {
  int n = resource_xptrs.size();
  LogicalVector was_using(n);
  for (int i = 0; i < n; ++i) {
    XPtr<DiscreteResource> ptr(VECTOR_ELT(resource_xptrs, i));
    was_using[i] = ptr->release_full(patient_id, amounts[i], purge_queues);
  }
  return was_using;
}

// [[Rcpp::export]]
LogicalVector discrete_resource_release_all_if_using_cpp(List resource_xptrs, int patient_id,
                                                          IntegerVector amounts) {
  int n = resource_xptrs.size();
  LogicalVector was_using(n);
  for (int i = 0; i < n; ++i) {
    XPtr<DiscreteResource> ptr(VECTOR_ELT(resource_xptrs, i));
    bool had = ptr->is_patient_using(patient_id);
    if (had) ptr->attempt_free_if_using(patient_id, false, amounts[i]);
    was_using[i] = had;
  }
  return was_using;
}

// ── Clone helper ─────────────────────────────────────────────────────────────

static SEXP get_dotptr_from_env(SEXP wrapper_env) {
  if (!Rf_isEnvironment(wrapper_env)) Rcpp::stop("Expected a resource wrapper environment.");
  Rcpp::Environment env(wrapper_env);
  if (!env.exists(".ptr")) Rcpp::stop("Wrapper is missing a valid '.ptr' external pointer.");
  SEXP v = env[".ptr"];
  if (TYPEOF(v) != EXTPTRSXP) Rcpp::stop("Wrapper is missing a valid '.ptr' external pointer.");
  return v;
}

// [[Rcpp::export]]
Rcpp::List discrete_resource_clone_xptrs_cpp(SEXP wrapper_env, int n = 1) {
  if (n <= 0) Rcpp::stop("n must be >= 1");

  SEXP xptr = get_dotptr_from_env(wrapper_env);
  Rcpp::XPtr<DiscreteResource> src(xptr);
  const int capacity       = src->size();
  const bool lifo          = src->get_is_lifo();
  const int max_queue      = src->get_max_queue_capacity();
  const bool allow_multi   = src->get_allow_multiple_queue();

  Rcpp::List out(n);
  try {
    for (int i = 0; i < n; ++i) {
      auto* p = new DiscreteResource(capacity, lifo, max_queue, allow_multi);
      out[i] = Rcpp::XPtr<DiscreteResource>(p, true);
    }
  } catch (const std::bad_alloc&) {
    Rcpp::stop("Memory allocation failed while cloning resources (n=%d).", n);
  } catch (const std::exception& e) {
    Rcpp::stop("Error while cloning resources: %s", e.what());
  } catch (...) {
    Rcpp::stop("Unknown error while cloning resources.");
  }
  return out;
}
