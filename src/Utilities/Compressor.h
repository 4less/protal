#include <iostream>
#include <string>
#include <filesystem>
#include <cstdlib>   // for std::system
#include <stdexcept> // for std::runtime_error

class Compressor {
public:
    // Compress a file in-place using pigz with the given number of threads
    static void compressInPlace(const std::filesystem::path& inputFile, int threads)
    {
        namespace fs = std::filesystem;

        // Validate inputs
        if (threads <= 0)
            throw std::invalid_argument("Thread count must be greater than zero.");

        if (inputFile.empty())
            throw std::invalid_argument("Input file path is empty.");

        if (!fs::exists(inputFile))
            throw std::invalid_argument("Input file does not exist: " + inputFile.string());

        // Construct the pigz command
        std::string command = "pigz -f -p " + std::to_string(threads)
                            + " " + inputFile.string();


        // Execute and wait for completion
        int ret = std::system(command.c_str());
        if (ret == -1) {
            throw std::runtime_error("Failed to invoke system() to run pigz.");
        } else if (ret != 0) {
            throw std::runtime_error("pigz exited with non-zero status: " + std::to_string(ret));
        }

        // Verify compression result
        fs::path gzFile = inputFile;
        gzFile += ".gz";

        if (!fs::exists(gzFile)) {
            throw std::runtime_error("pigz completed but output file not found: " + gzFile.string());
        }
    }
};