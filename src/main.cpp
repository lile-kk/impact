#include <iostream>
#include <string>

// #include "src/run/Smack.h"

int main(int argc, char** argv) {
    bool* arr = new bool[10];
    std::cout << arr << std::endl;
    return 0;
} 
/**
//test
void test_impact() {
    // Create an instance of Smack

	Smack smackinstance = new Smack();

	// Check if something is written after the command
	if (args.length > 0) {

		// Read input file and parse
		try {
			System.out.println("=========== Welcome to Impact =========="); 	
			System.out.println("= Impact is Alpha code. This means the =");
			System.out.println("= results are not validated and should =");
			System.out.println("= not be trusted for use.              ="); 
			System.out.println("=                                      =");
			System.out.println("= Impact is free software and comes    =");
			System.out.println("= with no warranties whatsoever        =");
			System.out.println("=                                      =");
			System.out.println("= Report any bugs at our webpage:      =");
			System.out.println("=============== ENJOY!!!! ==============\n\n");
			System.out.println("Processing file: " + args[0]);
			System.out.println("*** Initializing ***");
			smackinstance.initialize(args[0]);
			// Assemble the mass matrix
			System.out.println("*** Assembling the Mass Matrix ***");
			smackinstance.assembleMassMatrix();
			// Set up the prerequisites
			System.out.println("*** Setting Initial Conditions ***");
			smackinstance.setInitialConditions();
			// Solve the problem
			System.out.println("*** Initiating Solver ***");
			smackinstance.solve();
			// Post the results and clean up
			smackinstance.post();
		} catch (Error e) {
			System.out.println(e);
			return;
		}
	} else {
		System.out.println("Impact - A simple explicit dynamic program\n");
		System.out.println("Syntax: java run.Impact sourcefile.in");
	}
}

*/

