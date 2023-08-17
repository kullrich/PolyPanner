#ifndef RESOURCE_H
#define RESOURCE_H

#include <thread>
#include <condition_variable>
#include <chrono>
#include <atomic>

using namespace std;

class IntResourceManager {
private:
  // tracks current usage of resource
  int m_resource;
  
  condition_variable m_cv;
  mutex m_mtx;
  
public:
  IntResourceManager(int total_resource): m_resource(total_resource) {};

  // total resource
  void set_total_resource(int usage);
  int get_total_resource();

  // wait until resource is free
  void obtain_resource(int idx, int usage);
  void free_resource(int idx, int usage);
  
};

#endif
