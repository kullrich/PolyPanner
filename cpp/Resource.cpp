#include <iostream>
#include "Resource.h"

void IntResourceManager::obtain_resource(int idx, int usage)
{
    std::unique_lock<std::mutex> lk(m_mtx);

    // cout << "PRE_WAIT: " << idx << ", current_usage=" << get_total_resource() << ", req_usage=" << usage << endl;
    while (usage > m_resource)
      m_cv.wait(lk);
    // cout << "START: " << idx << ", current_usage=" << get_total_resource() << ", req_usage=" << usage << endl;
    m_resource -= usage;
}

void IntResourceManager::free_resource(int idx, int usage)
{
    std::unique_lock<std::mutex> lk(m_mtx);
    m_resource += usage;
    // cout << "DONE: " << idx << ", current_usage=" << get_total_resource() << ", req_usage=" << usage << endl;
    lk.unlock();
    m_cv.notify_one();
}

void IntResourceManager::set_total_resource(int resource)
{
  m_resource = resource;
}

int IntResourceManager::get_total_resource()
{
    return m_resource;
}
