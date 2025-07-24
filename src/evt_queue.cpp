// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
#include <queue>
#include <unordered_map>
#include <vector>
#include <string>
#include <stdexcept>
using namespace Rcpp;

struct Event {
  int patient_id;
  std::string event_name;
  double time;
  int priority; // lower is higher priority
  size_t version; // for tracking stale events
  
  Event(int pid, const std::string& name, double t, int p, size_t v) 
    : patient_id(pid), event_name(name), time(t), priority(p), version(v) {}
};

struct EventCompare {
  bool operator()(const Event& a, const Event& b) const {
    if (a.time != b.time) return a.time > b.time; // earlier time = higher priority
    if (a.priority != b.priority) return a.priority > b.priority; // lower index = higher priority
    return a.version > b.version; // deterministic ordering for tie-breaking
  }
};

struct Key {
  int patient_id;
  std::string event_name;
  
  bool operator==(const Key& other) const {
    return patient_id == other.patient_id && event_name == other.event_name;
  }
};

struct KeyHash {
  std::size_t operator()(const Key& k) const {
    // Improved hash function with better mixing
    std::size_t h1 = std::hash<int>()(k.patient_id);
    std::size_t h2 = std::hash<std::string>()(k.event_name);
    return h1 ^ (h2 << 1);
  }
};

struct EventInfo {
  double time;
  size_t version;
  
  // Add default constructor
  EventInfo() : time(0.0), version(0) {}
  
  // Keep your existing constructor
  EventInfo(double t, size_t v) : time(t), version(v) {}
};

class EventQueue {
private:
  std::priority_queue<Event, std::vector<Event>, EventCompare> heap;
  std::unordered_map<Key, EventInfo, KeyHash> lookup;
  std::unordered_map<std::string, int> event_priority;
  size_t next_version;
  size_t stale_count;
  static const size_t REBUILD_THRESHOLD = 10000;
  static const size_t MAX_STALE_ITEMS = 1000000; // ~58MB of stale events
  
  void clean_top() {
    while (!heap.empty()) {
      const Event& top_event = heap.top();
      Key k = {top_event.patient_id, top_event.event_name};
      auto it = lookup.find(k);
      
      if (it != lookup.end() && 
          it->second.time == top_event.time && 
          it->second.version == top_event.version) {
        break; // Found valid top event
      }
      
      heap.pop();
      if (stale_count > 0) stale_count--; // Prevent underflow
    }
  }
  
  void rebuild_if_needed() {
    if ((stale_count > REBUILD_THRESHOLD && stale_count > heap.size() / 2) ||
        stale_count > MAX_STALE_ITEMS) {
      rebuild_heap();
      
    }
  }
  
  void rebuild_heap() {
    std::vector<Event> valid_events;
    valid_events.reserve(lookup.size());
    
    for (const auto& pair : lookup) {
      const Key& k = pair.first;
      const EventInfo& info = pair.second;
      
      auto priority_it = event_priority.find(k.event_name);
      if (priority_it != event_priority.end()) {
        valid_events.emplace_back(k.patient_id, k.event_name, info.time, 
                                  priority_it->second, info.version);
      } else {
        // Defensive check: event in lookup but not in event_priority
        throw std::runtime_error("Inconsistent state: event '" + k.event_name + 
                                 "' found in lookup but not in priority table");
      }
    }
    
    // Replace heap with valid events only
    heap = std::priority_queue<Event, std::vector<Event>, EventCompare>(
      EventCompare(), std::move(valid_events));
    stale_count = 0;
  }
  
public:
  EventQueue(const std::vector<std::string>& priority_order) 
    : next_version(0), stale_count(0) {
    event_priority.reserve(priority_order.size());
    for (size_t i = 0; i < priority_order.size(); ++i) {
      event_priority[priority_order[i]] = static_cast<int>(i);
    }
  }
  
  void push(int patient_id, const std::string& event_name, double time) {
    auto priority_it = event_priority.find(event_name);
    if (priority_it == event_priority.end()) {
      throw std::runtime_error("Unknown event type: " + event_name);
    }
    
    Key k = {patient_id, event_name};
    size_t version = next_version++;
    int priority = priority_it->second;
    
    Event e(patient_id, event_name, time, priority, version);
    heap.push(e);
    
    // Update lookup, marking old version as stale if it exists
    auto lookup_it = lookup.find(k);
    if (lookup_it != lookup.end()) {
      stale_count++; // Old version becomes stale
    }
    lookup[k] = EventInfo(time, version);
    
    rebuild_if_needed();
  }
  
  void push_multiple(int patient_id, const std::vector<std::string>& event_names, 
                     const std::vector<double>& times) {
    if (event_names.size() != times.size()) {
      throw std::runtime_error("Event names and times vectors must have same length");
    }
    
    for (size_t i = 0; i < event_names.size(); ++i) {
      push(patient_id, event_names[i], times[i]);
    }
  }
  
  Event top() {
    clean_top();
    if (heap.empty()) {
      throw std::runtime_error("Queue is empty");
    }
    return heap.top();
  }
  
  void pop() {
    clean_top();
    if (heap.empty()) {
      throw std::runtime_error("Queue is empty");
    }
    
    const Event& e = heap.top(); // Store reference before popping
    Key k = {e.patient_id, e.event_name};
    lookup.erase(k);
    heap.pop();
  }
  
  Event pop_and_return() {
    Event e = top(); // This calls clean_top() and checks if empty
    Key k = {e.patient_id, e.event_name};
    lookup.erase(k);
    heap.pop();
    return e;
  }
  
  void remove(int patient_id, const std::string& event_name) {
    Key k = {patient_id, event_name};
    auto it = lookup.find(k);
    if (it != lookup.end()) {
      lookup.erase(it);
      stale_count++; // Event becomes stale in heap
      rebuild_if_needed();
    }
  }
  
  void remove_multiple(int patient_id, const std::vector<std::string>& event_names) {
    for (const auto& event_name : event_names) {
      remove(patient_id, event_name);
    }
  }
  
  void modify(int patient_id, const std::string& event_name, double new_time, 
              bool create_if_missing = false) {
    Key k = {patient_id, event_name};
    auto it = lookup.find(k);
    
    if (it != lookup.end()) {
      remove(patient_id, event_name);
      push(patient_id, event_name, new_time);
    } else if (create_if_missing) {
      push(patient_id, event_name, new_time);
    }
  }
  
  void modify_multiple(int patient_id, const std::vector<std::string>& event_names, 
                       const std::vector<double>& new_times, bool create_if_missing = false) {
    if (event_names.size() != new_times.size()) {
      throw std::runtime_error("Event names and times vectors must have same length");
    }
    
    for (size_t i = 0; i < event_names.size(); ++i) {
      modify(patient_id, event_names[i], new_times[i], create_if_missing);
    }
  }
  
  bool empty() const {
    return lookup.empty();
  }
  
  size_t size() const {
    return lookup.size();
  }
  
  bool has_event(int patient_id, const std::string& event_name) const {
    Key k = {patient_id, event_name};
    return lookup.find(k) != lookup.end();
  }
  
  std::vector<Event> peek_next(size_t n) {
    std::vector<Event> result;
    std::vector<Event> temp_storage;
    
    size_t actual_n = std::min(n, lookup.size());
    result.reserve(actual_n);
    temp_storage.reserve(actual_n);
    
    // Extract n events temporarily
    for (size_t i = 0; i < actual_n; ++i) {
      clean_top();
      if (heap.empty()) break;
      
      Event e = heap.top();
      result.push_back(e);
      temp_storage.push_back(e);
      
      // Remove from both heap and lookup temporarily
      Key k = {e.patient_id, e.event_name};
      lookup.erase(k);
      heap.pop();
    }
    
    // Restore all events with their original versions
    for (const auto& e : temp_storage) {
      Key k = {e.patient_id, e.event_name};
      lookup[k] = EventInfo(e.time, e.version);
      heap.push(e);
    }
    
    return result;
  }
  
  double get_event_time(int patient_id, const std::string& event_name) const {
    Key k = {patient_id, event_name};
    auto it = lookup.find(k);
    if (it == lookup.end()) {
      throw std::runtime_error("Event not found for this patient");
    }
    return it->second.time;
  }
  
  const std::unordered_map<Key, EventInfo, KeyHash>& get_lookup() const {
    return lookup;
  }
  
  const std::unordered_map<std::string, int>& get_event_priority() const {
    return event_priority;
  }
};

// Helper function to convert named vector to vectors of names and values
std::pair<std::vector<std::string>, std::vector<double>> 
  parse_named_vector(const NumericVector& named_vec) {
    std::vector<std::string> names;
    std::vector<double> values;
    
    CharacterVector vec_names = named_vec.names();
    if (vec_names.size() != named_vec.size() || Rf_isNull(vec_names)) {
      throw std::runtime_error("All elements in named vector must have names");
    }
    
    names.reserve(named_vec.size());
    values.reserve(named_vec.size());
    
    for (int i = 0; i < named_vec.size(); ++i) {
      names.push_back(as<std::string>(vec_names[i]));
      values.push_back(named_vec[i]);
    }
    
    return std::make_pair(names, values);
  }

// Expose via XPtr
typedef Rcpp::XPtr<EventQueue> EventQueuePtr;

// [[Rcpp::export]]
SEXP queue_create_cpp(std::vector<std::string> priority_order) {
  try {
    EventQueue* q = new EventQueue(priority_order);
    EventQueuePtr ptr(q, true);
    return ptr;
  } catch (const std::exception& e) {
    throw Rcpp::exception(e.what());
  }
}

// [[Rcpp::export]]
void new_event_cpp(SEXP ptr, int patient_id, NumericVector events) {
  try {
    EventQueuePtr q(ptr);
    auto parsed = parse_named_vector(events);
    q->push_multiple(patient_id, parsed.first, parsed.second);
  } catch (const std::exception& e) {
    throw Rcpp::exception(e.what());
  }
}

// [[Rcpp::export]]
List next_event_cpp(SEXP ptr, int n = 1) {
  try {
    EventQueuePtr q(ptr);
    
    if (n <= 0) {
      return List::create(
        _["patient_id"] = IntegerVector(),
        _["event_name"] = CharacterVector(),
        _["time"] = NumericVector()
      );
    }
    
    std::vector<Event> events = q->peek_next(static_cast<size_t>(n));
    
    std::vector<int> patient_ids;
    std::vector<std::string> event_names;
    std::vector<double> times;
    
    patient_ids.reserve(events.size());
    event_names.reserve(events.size());
    times.reserve(events.size());
    
    for (const auto& e : events) {
      patient_ids.push_back(e.patient_id);
      event_names.push_back(e.event_name);
      times.push_back(e.time);
    }
    
    return List::create(
      _["patient_id"] = patient_ids,
      _["event_name"] = event_names,
      _["time"] = times
    );
  } catch (const std::exception& e) {
    throw Rcpp::exception(e.what());
  }
}

// [[Rcpp::export]]
List next_event_pt_cpp(SEXP ptr, int patient_id, int n = 1) {
  try {
    EventQueuePtr q(ptr);
    const auto& lookup = q->get_lookup();               // Use getter
    const auto& event_priority = q->get_event_priority();  // Use getter
    
    if (n <= 0) {
      return List::create(
        _["patient_id"] = IntegerVector(),
        _["event_name"] = CharacterVector(),
        _["time"] = NumericVector()
      );
    }
    
    // Gather all events for the patient
    std::vector<Event> patient_events;
    for (const auto& pair : lookup) {
      const Key& k = pair.first;
      const EventInfo& info = pair.second;
      if (k.patient_id == patient_id) {
        auto prio_it = event_priority.find(k.event_name);
        int prio = (prio_it != event_priority.end()) ? prio_it->second : 0;
        patient_events.emplace_back(k.patient_id, k.event_name, info.time, prio, info.version);
      }
    }
    
    // If no events found, return empty
    if (patient_events.empty()) {
      return List::create(
        _["patient_id"] = IntegerVector(),
        _["event_name"] = CharacterVector(),
        _["time"] = NumericVector()
      );
    }
    
    // Sort events by time ascending, then priority ascending, then version ascending
    std::sort(patient_events.begin(), patient_events.end(), [](const Event& a, const Event& b) {
      if (a.time != b.time) return a.time < b.time;
      if (a.priority != b.priority) return a.priority < b.priority;
      return a.version < b.version;
    });
    
    size_t to_return = std::min(static_cast<size_t>(n), patient_events.size());
    
    // Prepare output vectors
    IntegerVector patient_ids(to_return);
    CharacterVector event_names(to_return);
    NumericVector times(to_return);
    
    for (size_t i = 0; i < to_return; ++i) {
      patient_ids[i] = patient_events[i].patient_id;
      event_names[i] = patient_events[i].event_name;
      times[i] = patient_events[i].time;
    }
    
    return List::create(
      _["patient_id"] = patient_ids,
      _["event_name"] = event_names,
      _["time"] = times
    );
  } catch (const std::exception& e) {
    throw Rcpp::exception(e.what());
  }
}


// [[Rcpp::export]]
void pop_event_cpp(SEXP ptr) {
  try {
    EventQueuePtr q(ptr);
    q->pop();
  } catch (const std::exception& e) {
    throw Rcpp::exception(e.what());
  }
}

// [[Rcpp::export]]
List pop_and_return_event_cpp(SEXP ptr) {
  try {
    EventQueuePtr q(ptr);
    Event e = q->pop_and_return();
    return List::create(
      _["patient_id"] = e.patient_id,
      _["event_name"] = e.event_name,
      _["time"] = e.time
    );
  } catch (const std::exception& e) {
    throw Rcpp::exception(e.what());
  }
}

// [[Rcpp::export]]
void remove_event_cpp(SEXP ptr, int patient_id, SEXP events) {
  try {
    EventQueuePtr q(ptr);
    
    if (TYPEOF(events) == STRSXP) {
      // Character vector
      CharacterVector event_names = as<CharacterVector>(events);
      std::vector<std::string> names;
      names.reserve(event_names.size());
      for (int i = 0; i < event_names.size(); ++i) {
        names.push_back(as<std::string>(event_names[i]));
      }
      q->remove_multiple(patient_id, names);
    } else if (TYPEOF(events) == REALSXP) {
      // Named numeric vector - check if it has names
      NumericVector named_vec = as<NumericVector>(events);
      if (Rf_isNull(named_vec.names())) {
        throw std::runtime_error("Numeric vector must have names when used for event removal");
      }
      auto parsed = parse_named_vector(named_vec);
      q->remove_multiple(patient_id, parsed.first);
    } else {
      throw std::runtime_error("Events must be character vector or named numeric vector");
    }
  } catch (const std::exception& e) {
    throw Rcpp::exception(e.what());
  }
}

// [[Rcpp::export]]
void modify_event_cpp(SEXP ptr, int patient_id, NumericVector events, 
                  bool create_if_missing = false) {
  try {
    EventQueuePtr q(ptr);
    auto parsed = parse_named_vector(events);
    q->modify_multiple(patient_id, parsed.first, parsed.second, create_if_missing);
  } catch (const std::exception& e) {
    throw Rcpp::exception(e.what());
  }
}

// [[Rcpp::export]]
bool queue_empty_cpp(SEXP ptr) {
  try {
    EventQueuePtr q(ptr);
    return q->empty();
  } catch (const std::exception& e) {
    throw Rcpp::exception(e.what());
  }
}

// [[Rcpp::export]]
int queue_size_cpp(SEXP ptr) {
  try {
    EventQueuePtr q(ptr);
    return static_cast<int>(q->size());
  } catch (const std::exception& e) {
    throw Rcpp::exception(e.what());
  }
}

// [[Rcpp::export]]
bool has_event_cpp(SEXP ptr, int patient_id, std::string event_name) {
  try {
    EventQueuePtr q(ptr);
    return q->has_event(patient_id, event_name);
  } catch (const std::exception& e) {
    throw Rcpp::exception(e.what());
  }
}


// [[Rcpp::export]]
double get_event_cpp(SEXP ptr, int patient_id, std::string event_name) {
  try {
    EventQueuePtr q(ptr);
    return q->get_event_time(patient_id, event_name);
  } catch (const std::exception& e) {
    throw Rcpp::exception(e.what());
  }
}