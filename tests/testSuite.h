#pragma once
#include <iostream>
#include <string>
#include <cmath>
#include <iomanip>

class TestSuite {
private:
    int totalTests_ = 0;
    int passedTests_ = 0;
    int failedTests_ = 0;
    std::string currentSection_;

public:
    void startSection(const std::string& sectionName) {
        currentSection_ = sectionName;
        std::cout << "\n========================================\n";
        std::cout << " RUNNING: " << sectionName << "\n";
        std::cout << "========================================\n";
    }

    void runTest(const std::string& testName, bool (*testFunc)()) {
        totalTests_++;
        std::cout << "  [RUNNING] " << testName << "... " << std::flush;
        if (testFunc()) {
            passedTests_++;
            std::cout << "\033[32m[PASS]\033[0m\n";
        } else {
            failedTests_++;
            std::cout << "\033[31m[FAIL]\033[0m\n";
        }
    }

    int printSummary() const {
        std::cout << "\n========================================\n";
        std::cout << "          TEST SUITE SUMMARY            \n";
        std::cout << "========================================\n";
        std::cout << " Total Suites Run: " << totalTests_ << "\n";
        std::cout << " Passed:          \033[32m" << passedTests_ << "\033[0m\n";
        std::cout << " Failed:          " << (failedTests_ > 0 ? "\033[31m" : "") 
                  << failedTests_ << (failedTests_ > 0 ? "\033[0m" : "") << "\n";
        std::cout << "----------------------------------------\n";
        
        if (failedTests_ == 0) {
            std::cout << " \033[32mALL TESTS PASSED SUCCESSFULLY!\033[0m\n";
        } else {
            std::cout << " \033[31mSOME TESTS FAILED!\033[0m\n";
        }
        std::cout << "========================================\n\n";

        return (failedTests_ == 0) ? 0 : 1;
    }
};

// Common assertion macros
#define TEST_ASSERT(cond, msg) \
    if (!(cond)) { \
        std::cerr << "\n    \033[31mAssertion Failed:\033[0m " << msg << " (Line " << __LINE__ << ")\n   "; \
        return false; \
    }

#define TEST_ASSERT_NEAR(val1, val2, eps, msg) \
    if (std::abs((val1) - (val2)) > (eps)) { \
        std::cerr << "\n    \033[31mAssertion Failed:\033[0m " << msg \
                  << " | Expected: " << (val2) << ", Got: " << (val1) << " (Line " << __LINE__ << ")\n   "; \
        return false; \
    }