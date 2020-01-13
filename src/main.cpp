#include <iostream>
#include <fstream>
#include <ctime>
#include <cstdio>
#include <csignal>
#include <getopt.h>
#include <sys/stat.h>

// New subcommand system provides all the subcommands that used to live here
#include "subcommand/subcommand.hpp"

using namespace std;
using namespace gyeet;

void gyeet_help(char** argv) {
    cerr << "gyeet: ultrafast sequence/graph aligner" << endl
         << endl
         << "usage: " << argv[0] << " <command> [options]" << endl
         << endl
         << gyeet::subcommand::PIPELINE << ":" << endl;
         
    gyeet::subcommand::Subcommand::for_each(gyeet::subcommand::PIPELINE, [](const gyeet::subcommand::Subcommand& command) {
        // Announce every subcommand we have
        
        // Pad all the names so the descriptions line up
        string name = command.get_name();
        name.resize(14, ' ');
        cerr << "  -- " << name << command.get_description() << endl;
     });
     
     cerr << endl << "For more commands, type `gyeet help`." << endl;
 }

// We make sure to compile main for the lowest common denominator architecture.
// This works on GCC and Clang. But we have to decalre main and then define it.
int main(int argc, char *argv[]) __attribute__((__target__("arch=x86-64")));

int main(int argc, char *argv[]) {

    // set a higher value for tcmalloc warnings
    setenv("TCMALLOC_LARGE_ALLOC_REPORT_THRESHOLD", "1000000000000000", 1);

    if (argc == 1) {
        gyeet_help(argv);
        return 1;
    }
    
    auto* subcommand = gyeet::subcommand::Subcommand::get(argc, argv);
    if (subcommand != nullptr) {
        // We found a matching subcommand, so run it
        return (*subcommand)(argc, argv);
    } else {
        // No subcommand found
        string command = argv[1];
        cerr << "error:[dg] command " << command << " not found" << endl;
        gyeet_help(argv);
        return 1;
    }

}
