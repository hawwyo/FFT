 #ifndef __TIME_MEASURE__HPP_INCLUDED__
 #define __TIME_MEASURE__HPP_INCLUDED__

 #include <ostream>
 #include <chrono>

 class time_measure {
    std::chrono::time_point<std::chrono::steady_clock> start, finish;
 public:
    time_measure() : start(std::chrono::steady_clock::now()) {
    }
    void stop() {
        finish = std::chrono::steady_clock::now();
    }
    int get() const {
    using namespace std::chrono;
        return duration_cast<milliseconds>(finish - start).count();
    }
 };

 std::ostream& operator<<(std::ostream& out, const time_measure& tm) {
    out << tm.get() << " ms";
    return out;
 }

 #endif
