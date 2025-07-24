#include <Rcpp.h>
#include <queue>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>

using namespace Rcpp;

// Structure to hold patient information in queue
struct QueuedPatient {
  int patient_id;
  int priority;
  int insertion_order;
  double queue_start_time;
  bool is_valid;  // For lazy deletion
  
  QueuedPatient(int id, int prio, int order, double start_time) : 
    patient_id(id), priority(prio), insertion_order(order), 
    queue_start_time(start_time), is_valid(true) {}
};

// Comparator for priority queue (higher priority first, then FIFO for same priority)
struct QueueComparator {
  bool operator()(const QueuedPatient& a, const QueuedPatient& b) const {
    if (a.priority != b.priority) {
      return a.priority < b.priority; // Higher priority first
    }
    return a.insertion_order > b.insertion_order; // FIFO for same priority
  }
};

// Structure to hold patient information when using resource
struct UsingPatient {
  int patient_id;
  double start_time;
  
  UsingPatient(int id, double time) : patient_id(id), start_time(time) {}
};

class DiscreteResource {
private:
  int total_capacity;
  int next_insertion_order;
  int current_max_priority;
  int total_entries_ever_added;
  int current_valid_entries;
  int operations_since_cleanup;
  
  std::vector<UsingPatient> patients_using;
  std::priority_queue<QueuedPatient, std::vector<QueuedPatient>, QueueComparator> patient_queue;
  
  std::unordered_map<int, int> patient_queue_count;
  std::unordered_map<int, std::vector<double>> patient_queue_start_times;
  
  static const int MAX_QUEUE_ENTRIES_PER_PATIENT = 1000;
  static const int CLEANUP_FREQUENCY = 10000;
  
  void cleanup_queue_top() {
    while (!patient_queue.empty() && 
           (!patient_queue.top().is_valid || 
           patient_queue_count.find(patient_queue.top().patient_id) == patient_queue_count.end() ||
           patient_queue_count[patient_queue.top().patient_id] == 0)) {
      patient_queue.pop();
    }
  }
  
  void full_queue_cleanup() {
    std::vector<QueuedPatient> valid_patients;
    
    while (!patient_queue.empty()) {
      QueuedPatient patient = patient_queue.top();
      patient_queue.pop();
      
      if (patient_queue_count.find(patient.patient_id) != patient_queue_count.end() &&
          patient_queue_count[patient.patient_id] > 0) {
        patient.is_valid = true;
        valid_patients.push_back(patient);
      }
    }
    
    for (const auto& patient : valid_patients) {
      patient_queue.push(patient);
    }
    
    total_entries_ever_added = static_cast<int>(valid_patients.size());
    current_valid_entries = static_cast<int>(valid_patients.size());
  }
  
  void check_and_cleanup() {
    operations_since_cleanup++;
    
    bool periodic_cleanup = (operations_since_cleanup >= CLEANUP_FREQUENCY);
    bool threshold_cleanup = (total_entries_ever_added - current_valid_entries > 
                                static_cast<int>(current_valid_entries * 0.5));
    
    if (periodic_cleanup || threshold_cleanup) {
      full_queue_cleanup();
      operations_since_cleanup = 0;
    }
  }
  
public:
  DiscreteResource(int n) : 
  total_capacity(n), 
  next_insertion_order(0),
  current_max_priority(1),
  total_entries_ever_added(0),
  current_valid_entries(0),
  operations_since_cleanup(0) {
    if (n < 0) {
      stop("Resource capacity must be >= 0");
    }
    if (n > 0) {
      patients_using.reserve(n);
    }
  }
  
  int size() const { return total_capacity; }
  int queue_size() const { return current_valid_entries; }
  int n_free() const { return total_capacity - patients_using.size(); }
  
  std::vector<int> patients_using_ids() const {
    std::vector<int> ids;
    ids.reserve(patients_using.size());
    for (const auto& patient : patients_using) {
      ids.push_back(patient.patient_id);
    }
    return ids;
  }
  
  std::vector<double> patients_using_times() const {
    std::vector<double> times;
    times.reserve(patients_using.size());
    for (const auto& patient : patients_using) {
      times.push_back(patient.start_time);
    }
    return times;
  }
  
  bool is_patient_in_queue(int patient_id) const {
    return patient_queue_count.find(patient_id) != patient_queue_count.end() &&
      patient_queue_count.at(patient_id) > 0;
  }
  
  bool is_patient_using(int patient_id) const {
    return std::find_if(patients_using.begin(), patients_using.end(),
                        [patient_id](const UsingPatient& p) {
                          return p.patient_id == patient_id;
                        }) != patients_using.end();
  }
  
  bool attempt_block(int patient_id, int priority, double start_time) {
    if (patient_queue_count.find(patient_id) != patient_queue_count.end() &&
        patient_queue_count[patient_id] >= MAX_QUEUE_ENTRIES_PER_PATIENT) {
      stop("Patient exceeds maximum queue entries limit");
    }
    
    if (current_valid_entries == 0 && n_free() > 0) {
      patients_using.emplace_back(patient_id, start_time);
      return true;
    }
    
    if (n_free() > 0 && current_valid_entries > 0) {
      cleanup_queue_top();
      
      if (!patient_queue.empty()) {
        QueuedPatient next_in_line = patient_queue.top();
        if (next_in_line.patient_id == patient_id && 
            patient_queue_count.find(patient_id) != patient_queue_count.end() &&
            patient_queue_count[patient_id] > 0) {
          
          patient_queue.pop();
          patients_using.emplace_back(patient_id, start_time);
          
          patient_queue_count[patient_id]--;
          if (patient_queue_count[patient_id] == 0) {
            patient_queue_count.erase(patient_id);
            patient_queue_start_times.erase(patient_id);
          } else {
            patient_queue_start_times[patient_id].erase(patient_queue_start_times[patient_id].begin());
          }
          current_valid_entries--;
          
          return true;
        }
      }
    }
    
    patient_queue.emplace(patient_id, priority, next_insertion_order++, start_time);
    
    if (patient_queue_count.find(patient_id) == patient_queue_count.end()) {
      patient_queue_count[patient_id] = 0;
      patient_queue_start_times[patient_id].reserve(10);
    }
    
    patient_queue_count[patient_id]++;
    patient_queue_start_times[patient_id].push_back(start_time);
    
    current_max_priority = std::max(current_max_priority, priority);
    
    total_entries_ever_added++;
    current_valid_entries++;
    
    check_and_cleanup();
    
    return false;
  }
  
  void attempt_free(int patient_id, bool remove_all = false) {
    bool found_in_using = false;
    
    if (remove_all) {
      auto original_size = patients_using.size();
      patients_using.erase(
        std::remove_if(patients_using.begin(), patients_using.end(),
                       [patient_id](const UsingPatient& p) {
                         return p.patient_id == patient_id;
                       }),
                       patients_using.end());
      found_in_using = original_size != patients_using.size();
    } else {
      auto it = std::find_if(patients_using.begin(), patients_using.end(),
                             [patient_id](const UsingPatient& p) {
                               return p.patient_id == patient_id;
                             });
      if (it != patients_using.end()) {
        patients_using.erase(it);
        found_in_using = true;
      }
    }
    
    if (found_in_using) {
      return;
    }
    
    if (patient_queue_count.find(patient_id) != patient_queue_count.end() &&
        patient_queue_count[patient_id] > 0) {
      
      patient_queue_count[patient_id]--;
      current_valid_entries--;
      
      if (!patient_queue_start_times[patient_id].empty()) {
        patient_queue_start_times[patient_id].erase(patient_queue_start_times[patient_id].begin());
      }
      
      if (patient_queue_count[patient_id] == 0) {
        patient_queue_count.erase(patient_id);
        patient_queue_start_times.erase(patient_id);
      }
      
      check_and_cleanup();
    }
  }
  
  void attempt_free_if_using(int patient_id, bool remove_all = false) {
    if (remove_all) {
      patients_using.erase(
        std::remove_if(patients_using.begin(), patients_using.end(),
                       [patient_id](const UsingPatient& p) {
                         return p.patient_id == patient_id;
                       }),
                       patients_using.end());
    } else {
      auto it = std::find_if(patients_using.begin(), patients_using.end(),
                             [patient_id](const UsingPatient& p) {
                               return p.patient_id == patient_id;
                             });
      if (it != patients_using.end()) {
        patients_using.erase(it);
      }
    }
  }
  
  std::vector<int> next_patient_in_line(int n = 1) {
    std::vector<int> result;
    std::priority_queue<QueuedPatient, std::vector<QueuedPatient>, QueueComparator> temp_queue;
    
    int items_to_check = std::min(n * 2, static_cast<int>(patient_queue.size()));
    std::priority_queue<QueuedPatient, std::vector<QueuedPatient>, QueueComparator> original_queue = patient_queue;
    
    for (int i = 0; i < items_to_check && !original_queue.empty(); ++i) {
      temp_queue.push(original_queue.top());
      original_queue.pop();
    }
    
    int count = std::min(n, current_valid_entries);
    while (static_cast<int>(result.size()) < count && !temp_queue.empty()) {
      QueuedPatient patient = temp_queue.top();
      temp_queue.pop();
      
      if (patient_queue_count.find(patient.patient_id) != patient_queue_count.end() &&
          patient_queue_count.at(patient.patient_id) > 0) {
        result.push_back(patient.patient_id);
      }
    }
    
    return result;
  }
  
  std::vector<double> queue_start_times() {
    std::vector<double> result;
    result.reserve(current_valid_entries);
    
    std::vector<int> patient_ids = next_patient_in_line(current_valid_entries);
    
    std::unordered_map<int, int> patient_index;
    
    for (int patient_id : patient_ids) {
      int& index = patient_index[patient_id];
      
      if (patient_queue_start_times.find(patient_id) != patient_queue_start_times.end() &&
          index < static_cast<int>(patient_queue_start_times.at(patient_id).size())) {
        result.push_back(patient_queue_start_times.at(patient_id)[index]);
      } else {
        result.push_back(0.0);
      }
      index++;
    }
    
    return result;
  }
  
  void modify_priority(int patient_id, int new_priority) {
    if (patient_queue_count.find(patient_id) == patient_queue_count.end() ||
        patient_queue_count[patient_id] == 0) {
      return;
    }
    
    int count = patient_queue_count[patient_id];
    std::vector<double> start_times = patient_queue_start_times[patient_id];
    
    for (int i = 0; i < count; ++i) {
      patient_queue.emplace(patient_id, new_priority, next_insertion_order++, start_times[i]);
      total_entries_ever_added++;
    }
    
    check_and_cleanup();
  }
  
  void add_resource(int n_to_add) {
    if (n_to_add <= 0) {
      stop("n_to_add must be positive");
    }
    total_capacity += n_to_add;
  }
  
  void remove_resource(int n_to_remove, double current_time) {
    if (n_to_remove <= 0) {
      stop("n_to_remove must be positive");
    }
    
    if (n_to_remove > total_capacity) {
      stop("Cannot remove more resources than available");
    }
    
    total_capacity -= n_to_remove;
    
    while (static_cast<int>(patients_using.size()) > total_capacity) {
      UsingPatient patient = patients_using.back();
      patients_using.pop_back();
      
      int new_priority = current_max_priority + 1;
      current_max_priority = new_priority;
      
      patient_queue.emplace(patient.patient_id, new_priority, next_insertion_order++, current_time);
      
      if (patient_queue_count.find(patient.patient_id) == patient_queue_count.end()) {
        patient_queue_count[patient.patient_id] = 0;
        patient_queue_start_times[patient.patient_id].reserve(10);
      }
      
      patient_queue_count[patient.patient_id]++;
      patient_queue_start_times[patient.patient_id].push_back(current_time);
      
      total_entries_ever_added++;
      current_valid_entries++;
    }
    
    check_and_cleanup();
  }
};

// XPtr constructor function
// [[Rcpp::export]]
SEXP create_discrete_resource_cpp(int n) {
  DiscreteResource* ptr = new DiscreteResource(n);
  XPtr<DiscreteResource> xptr(ptr, true);
  return xptr;
}

void validate_xptr(SEXP xptr) {
  if (TYPEOF(xptr) != EXTPTRSXP) {
    stop("Invalid external pointer");
  }
}

// [[Rcpp::export]]
int discrete_resource_size_cpp(SEXP xptr) {
  validate_xptr(xptr);
  XPtr<DiscreteResource> ptr(xptr);
  return ptr->size();
}

// [[Rcpp::export]]
int discrete_resource_queue_size_cpp(SEXP xptr) {
  validate_xptr(xptr);
  XPtr<DiscreteResource> ptr(xptr);
  return ptr->queue_size();
}

// [[Rcpp::export]]
int discrete_resource_n_free_cpp(SEXP xptr) {
  validate_xptr(xptr);
  XPtr<DiscreteResource> ptr(xptr);
  return ptr->n_free();
}

// [[Rcpp::export]]
IntegerVector discrete_resource_patients_using_cpp(SEXP xptr) {
  validate_xptr(xptr);
  XPtr<DiscreteResource> ptr(xptr);
  std::vector<int> ids = ptr->patients_using_ids();
  return wrap(ids);
}

// [[Rcpp::export]]
NumericVector discrete_resource_patients_using_times_cpp(SEXP xptr) {
  validate_xptr(xptr);
  XPtr<DiscreteResource> ptr(xptr);
  std::vector<double> times = ptr->patients_using_times();
  return wrap(times);
}

// [[Rcpp::export]]
bool discrete_resource_is_patient_in_queue_cpp(SEXP xptr, int patient_id) {
  validate_xptr(xptr);
  XPtr<DiscreteResource> ptr(xptr);
  return ptr->is_patient_in_queue(patient_id);
}

// [[Rcpp::export]]
bool discrete_resource_is_patient_using_cpp(SEXP xptr, int patient_id) {
  validate_xptr(xptr);
  XPtr<DiscreteResource> ptr(xptr);
  return ptr->is_patient_using(patient_id);
}

// [[Rcpp::export]]
bool discrete_resource_attempt_block_cpp(SEXP xptr, int patient_id, int priority, double start_time) {
  validate_xptr(xptr);
  XPtr<DiscreteResource> ptr(xptr);
  return ptr->attempt_block(patient_id, priority, start_time);
}

// [[Rcpp::export]]
void discrete_resource_attempt_free_cpp(SEXP xptr, int patient_id, bool remove_all = false) {
  validate_xptr(xptr);
  XPtr<DiscreteResource> ptr(xptr);
  ptr->attempt_free(patient_id, remove_all);
}

// [[Rcpp::export]]
void discrete_resource_attempt_free_if_using_cpp(SEXP xptr, int patient_id, bool remove_all = false) {
  validate_xptr(xptr);
  XPtr<DiscreteResource> ptr(xptr);
  ptr->attempt_free_if_using(patient_id, remove_all);
}

// [[Rcpp::export]]
IntegerVector discrete_resource_next_patient_in_line_cpp(SEXP xptr, int n = 1) {
  validate_xptr(xptr);
  XPtr<DiscreteResource> ptr(xptr);
  std::vector<int> patients = ptr->next_patient_in_line(n);
  return wrap(patients);
}

// [[Rcpp::export]]
NumericVector discrete_resource_queue_start_times_cpp(SEXP xptr) {
  validate_xptr(xptr);
  XPtr<DiscreteResource> ptr(xptr);
  std::vector<double> times = ptr->queue_start_times();
  return wrap(times);
}

// [[Rcpp::export]]
void discrete_resource_modify_priority_cpp(SEXP xptr, int patient_id, int new_priority) {
  validate_xptr(xptr);
  XPtr<DiscreteResource> ptr(xptr);
  ptr->modify_priority(patient_id, new_priority);
}

// [[Rcpp::export]]
void discrete_resource_add_resource_cpp(SEXP xptr, int n_to_add) {
  validate_xptr(xptr);
  XPtr<DiscreteResource> ptr(xptr);
  ptr->add_resource(n_to_add);
}

// [[Rcpp::export]]
void discrete_resource_remove_resource_cpp(SEXP xptr, int n_to_remove, double current_time) {
  validate_xptr(xptr);
  XPtr<DiscreteResource> ptr(xptr);
  ptr->remove_resource(n_to_remove, current_time);
}