// COSC3500/7502 Matrix multiply assignment

//The header file that defines the marking rubric
#include <rubric.h>

//misc. standard C/C++ libraries
#include <cstdlib>

//The main function simply interprets the command-line arguments, and passes them to GradeBot
int main(int argc, char** argv)
{
    //Default matrix dimension
    int N = DEFAULT_N;
    int threadCount = DEFAULT_THREADCOUNT;

    //Defines whether the CPU, GPU, and MPI benchmarks should be run
    int isEnabled[3] = { 1,1,1 };

    //Reading in the command-line arguments
    //first argument will the the executable path/name itself
    int* argsIn = 0;
    argsIn = (int*)std::malloc(sizeof(int) * argc);
    for (int i=0;i<argc;i++)
    {
        argsIn[i] = 0;
    }

    //Read in the command-line arguments
    for (int argIdx = 1; argIdx < argc; argIdx++)
    {
        int v = std::atoi(argv[argIdx]);
        argsIn[argIdx - 1] = v;

        if (argIdx == 1 && v > 0)
        {
            N = v;
        }
        if (argIdx == 2)
        {
            threadCount = v;
        }

        if (argIdx == 3)
        {
            isEnabled[RUBRIC_CPU] = v;
        }
        if (argIdx == 4)
        {
            isEnabled[RUBRIC_GPU] = v;
        }
        if (argIdx == 5)
        {
            isEnabled[RUBRIC_MPI] = v;
        }
    }

    int* args = 0;
    int argCount = 0;

    if (argc > 6)
    {
        args = &argsIn[5];
        argCount = argc - 6;
    }

    //Run the GradeBot routine
    GradeBot(N, threadCount, &isEnabled[0],args,argCount);

    std::free(argsIn);

}
