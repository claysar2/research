'''
This is an annotated example of a python-based MPI program that can be run on multiple 
processors (referred to as MPI tasks here).  

It creates the same 2D x, y, and HI arrays on each MPI task, and initializes them to 
the same values: theta and phi values for x and y, and then zeros for HI.  Each MPI
task then fills in its own PART of the HI array, which is combined together at the end 
using the MPI "Reduce" method.

Note: you may need to install mpi4py (http://mpi4py.scipy.org/docs/, 
"pip install mpi4py" if you're using Anaconda) to be able to run this.  Look at this
Git repository of mpi4py examples as well:  https://github.com/jbornschein/mpi4py-examples

To run this script, you'd type "mpirun -np 4 python MPI_example.py", where you can replace the 
4 after "-np" with any number of MPI tasks -- but make sure that the array you create divides 
evenly into that number of tasks!  (i.e., if you have 64x64 array for the Aitoff projections,
make sure that your number of MPI tasks divides into 64*64=4096, so use a power of 2!
'''

import numpy as np

# imports mpi4py's primary module and starts up MPI
from mpi4py import MPI

'''
Creates the MPI "communicator" and starts up MPI.
Note that "comm" is a Python object that you can refer to 
and get lots of information out of, like the total number
of MPI tasks and the current task's ID number
'''
comm = MPI.COMM_WORLD

'''
comm.rank tells you what your 'rank', or ID number, is.
comm.rank = 0 is the "root" process, which most people treat 
as the boss process
'''
if comm.rank == 0:
    print("I'm the root task!", comm.rank)
else:
    print("I'm NOT the root task!", comm.rank)

'''
We're setting up the Numpy meshgrid here.  
This is being done on every processor, so it's duplicated 
on every MPI task.  We won't modify the values (though we will
reshape the array), so we don't actually do anything to communicate
it later.
Note: x is theta, y is phi in our all-sky maps.
'''
Dx, Dy = 0.1, 0.05
x, y = np.mgrid[slice(-np.pi, np.pi + Dx, Dx),
                slice(-np.pi/2, np.pi/2 + Dy, Dy) ]


print("Task", comm.rank, "x",x.shape,x.size)
print("Task", comm.rank, "y",y.shape,y.size)

'''
we're going to reshape this array into a 1D array so we can loop over
it more easily, so we want to keep track of its original shape
'''
original_shape = x.shape

# make the arrays 1D
x = np.reshape(x,-1)
y = np.reshape(y,-1)

# make an HI array that is the same size as our x,y arrays (to store data)
HI = np.zeros_like(x)


print("Task", comm.rank, "x new",x.shape,x.size)
print("Task", comm.rank, "y new",y.shape,y.size)

'''
Based on the size of the array and the number of
MPI tasks, figure out how many array elements each MPI task
is responsible for as we step through the array.  Note that 
the // explicitly means do integer arithmetic, so dN must be
an int.
'''
dN = x.size // comm.size

'''
When I calculate dN, I am implicitly assuming that the array size and number of
MPI tasks divide evenly.  If they do NOT do that, we'll get
weird results, so check.
'''
if (x.size % comm.size == 0):
    print("Task", comm.rank, "no leftovers!")
else:
    print("Task", comm.rank, "your array size and MPI tasks do not divide evenly!")
    comm.Abort(errorcode=123)

'''
Calculate start and end indices for loop on THIS TASK
based on the rank.  Since the ranks go from 0 through N-1,
we can do some fairly straightforward arithmetic to ensure
that each rank does dN tasks.
'''
start = comm.rank*dN
end = (comm.rank+1)*dN

# print out the indices for this task
print("I am on rank", comm.rank, "and my start and end indices are:", start, end, "of", x.size)

'''
In this loop, each MPI task will loop over only its part of the HI array, leaving the rest of 
the array set to zero.  This is really important, and analogous to what you're going to want 
to try doing with yt - each MPI task will open the dataset on its own (once!) and then loop
over its part of the arrays.
'''
for i in range(start,end):
    HI[i] = comm.rank + x[i] + y[i]

'''
We now have to recombine all of the partial HI arrays into a single array.  We're going to 
do this on MPI task 0, the "root" task.  So, we need an array that we're going to combine all
of the data into, which we will initialize here.  It MUST be the same size as the HI array!
'''
reduce_array = np.zeros_like(HI)

'''
MPI reduce operation: takes the partial HI arrays on all of
the individual MPI tasks, copies them to the root task (comm.rank = 0),
adds them upp, and stores them in "reduce_array".  Note that we have
been careful to make each MPI task loop over only its part of the HI array,
which was set to zero everwhere before we started, so we aren't losing or 
duplicating any data.  There are lots of operations you can do here, but summing
things is a pretty typical one.
'''
comm.Reduce(HI,reduce_array,op=MPI.SUM)

# print out just to see it did something (on the root task)
if comm.rank == 0:
    print("on rank 0 my array is",reduce_array)

'''
reshape all of the arrays to be our original 2D arrays, which is what we 
need for Aitoff projections.  We will reshape x and y in place, but we'll
put reduce_array back into the HI array.  Note that while this happens on
all of the MPI tasks, ONLY the root task (rank 0) has the complete HI array!
'''

x = np.reshape(x,original_shape)
y = np.reshape(y,original_shape)
HI = np.reshape(reduce_array,original_shape)

# print it out on the root task
if comm.rank == 0:
    print("On task zero:  x, y, HI shapes:", x.shape, y.shape, HI.shape)



