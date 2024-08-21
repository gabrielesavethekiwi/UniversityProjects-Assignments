# Quantum Physics Course Weekly Assignments

I am collecting posthumously the weekly assignments from the course "Quantum Information and Computing". 

It has been a pivotal step in my development as a scientist, as I had limited exposure to programming and computational techniques before this course.

Through this course, I was introduced to the following categorical imperative for scientific software:

**Scientific software must satisfy the following priorities in descending order:**
1. **Correctness**
2. **Numerical stability**
3. **Accurate discretization**
4. **Flexibility**
5. **Efficiency**
6. **(Portability)**

Let me briefly comment:
1. Correctness - Its importance is exacerbated by dealing with quantum mechanics. We can include an automatic proof of correctness in the code  through pre-conditions, which must be true before a method is called, post-conditions, which must be true after a method is called, and invariants, which must always be true. If we imagine our method interfacing with a client metaphorically, a pre-condition tells the client "this is what I expect from you", a post-condition means "this is what I promise to do for you" and an invariant corresponds to "if this was true before you called me (the method), I promise it will still be true once I'm done". 
2. **Numerical stability** - Be cautious when working with real numbers, as they are prone to approximation errors. For example, a condition like x == y is unreliable due to the inherent imprecision of real number representations. Additionally, be aware of issues like catastrophic cancellation and error propagation, particularly in iterative methods where errors can grow exponentially.
3. **Accurate discretization** - The discretization must respect the symmetries of the system being simulated. Furthermore, let \( T \) be an observable defined on a certain interval, and let \( T_h \) be its discretized version computed on a grid with spacing \( h \). Then, as \( h \) approaches zero, \( T_h \) must converge to \( T \):

\[
\lim_{h \to 0} T_h = T
\]

5. **Efficiency** - If you took care of the previous steps, it's now much easier to consider the leading order algorithmic complexity and improve the efficiency of your program. "95% of time is spent on 5% of the code"


## Course Overview

From the course webpage [Quantum Information and Computation](https://didattica.unipd.it/off/2020/LM/SC/SC2443/000ZZ/SCP8082721/N0):


> "The course aims to introduce the students to computational quantum physics,  
> one of the most versatile simulation approaches exploited in quantum science.  
> It will provide a hands-on introduction to these methods and will present some  
> of tensor network methods' most successful and promising applications.  
> Indeed, they are routinely used to characterize low-dimensional equilibrium and  
> out-of-equilibrium quantum processes to guide and support the development of  
> quantum technologies."

