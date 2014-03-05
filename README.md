<!-- File generated by tutorialj -->
# conifer

Installing from source
----------------------

Requires: gradle, git, eclipse

- Clone the repository
- Type ``gradle eclipse`` from the root of the repository
- From eclipse:
  - ``Import`` in ``File`` menu
  - ``Import existing projects into workspace``
  - Select the root of the newly created repo
  - Deselect ``Copy projects into workspace`` to avoid having duplicates



Bayesian inference using Blang
-----------------------------

You will first learn to use ``Blang``, a new package where you can
write Bayesian models in a way reminiscent to JAGS, BUGS, etc. The 
main advantage of Blang is that you can easily add other types of 
random variables. This will be important in this exercise as we will
need tree-valued random variables.

### Small example: predicting the number of future members of the human species

Before creating new types of random variables, let us start with an example using
familiar, real-valued random variable to familiarize ourself with the main Blang
concepts. To do so, we will follow the lines of the most basic version of the 
[Doomsday argument](http://en.wikipedia.org/wiki/Doomsday_argument).

#### Overview

Suppose you want to predict the number N of future members of the human species. 

From archaeological records, we know there are about 0.06 trillion humans have been 
born so far. Call this our observation Y. Let us assume that Y is uniformly 
distributed between 0 and N.

So if we put a prior on N, we can compute a posterior on the quantity of interest,
N - Y. Let us be optimistic, and let us put an exponential prior on N with a mean of
10.0 trillion humans.

#### Basic concepts in Blang

Let us look now at how we can write this model in Blang.


<sub>From:[conifer.BlangSmallExample](src/test/java//conifer/BlangSmallExample.java)</sub>

Blang provides a way of building a factor graph declaratively. 

Recall that factor graphs are bipartite undirected graphs with two types 
of nodes: **factors** and **variables**. Each factor computes one factor
in a target unnormalized density broken into a product of factors. 
Each variable holds one coordinate of a state. 

For now, think of a factor as a prior or conditional probability density 
(although Blang can be used in other contexts such as for Markov random fields, 
we will focus on factor graphs
induced from directed models in this tutorial).

**Factors** are declared 
using the ``@DefineFactor`` annotation (an annotation is a mechanism designed 
to extend the functionality of the Java language). Blang will 
look for fields with a DefineFactor annotations and build model composed of 
these factors. 

We will look at the inner working of factors later on,
for now we focus on how to use them. For example, in our example, there are two 
factors, one for the prior (an exponential density), and one for the likelihood
(a uniform density).
 
**Variables** are fields within factors classes that have a 
``@FactorArgument`` and ``@FactorComponent`` annotations (the differences between 
the two will be explained shortly). These fields can be of any type (more on that
later), and they typically correspond to either the parameters or the realization 
of a probability density.

In order to use existing factors, the only thing to do with variables is
to *connect* them across factors. In our example, this is achieved via
``Exponential.on(likelihood.parameters.max)``. Here ``Exponential.on(x)`` is just a 
a static function short-hand for creating a unit mean exponential density with 
realization x, 
``new Exponential(x, ...)``. The consequence of doing this is that the realization field
of the exponential prior is ``==`` (object equality, i.e. two references to the same
location in the heap) to the realization field of 
the maximum parameter of the uniform density factor. Blang automatically discovers 
these relationships and use them to build edges in the factor graph.

#### Building an MCMC sampler

Obtaining an MCMC sampler from a model as specified above is very simple:

1. Make the class in which you declared the factor extend ``MCMCRunner``.
2. In the main function, create an instance of that class, and call ``.run()``.

If you run the sampler, you will see an histogram created in ``posteriorOnMax.pdf``,
showing the approximate posterior over the number of future members of the human 
species under our simple model.

#### Differentiating parameters from realisations

As mentioned earlier, there are two types of factor arguments in our setup: 
parameters and realisations. As we will see shortly, it 
is useful to differentiate the two. This is done via the syntax 
``@FactorArgument(makeStochastic = true)``, which indicates that the argument
is a realisation (therefore making the variable stochastic/random).

See ``bayonet.distributions.Exponential.java`` and ``bayonet.distributions.Uniform.java``
for examples of this construction.

#### Specifying what should be sampled

All the latent (non-observed) random variables are sampled by default. Blang uses the 
following strategy to determine the random variables: the variable should appear as annotated by  
``@FactorArgument(makeStochastic = true)`` in one of the factors (see above).

A random variable can be either latent or observed. To mark a variable as observed,
simply use ``@DefineFactor(onObservations = true)`` on the factor in which the 
observed variable is declared as stochastic.

#### Processing the samples

The simplest way of process the samples is to override the function ``process()``
as follows:



```java
protected void process(blang.processing.ProcessorContext)
{
    samples.add(likelihood.parameters.max.getValue() - observation.getValue());
}
```
<sub>From:[conifer.BlangSmallExample](src/test/java//conifer/BlangSmallExample.java)</sub>

### Extending Blang

Let us see now Blang can be extended beyond its (currently very minimal) set of built-in 
factor and variables. At the same time, some more advanced programming concepts will be 
introduced (generics, test units, annotations). 

#### Creating a new factor

In this part, you will create a new type of factor to model Gamma densities.

Here are the main steps:

1. Create a class ``Gamma``. Make this class implement ``Factor``.
2. Create three fields of type ``RealVariable``, and annotated with 
   ``@FactorArgument``: one for the realization, one for the shape, one for the rate.
   Note that these should also be ``public`` and ``final``. The keyword ``public`` simply 
   means they are accessible from functions outside of the current one. 
   The keyword ``final`` means that the value of the field cannot be changed.
3. Implement ``logDensity()``. You can the existing implementation provided by 
   the math-commons project (see how I did it in the ``Exponential`` class).
4. Create a constructor.

That's it! We now have a working factor.

Before we move on, just a quick not on why we made the fields ``final``. Recall that almost all 
variables in java are in fact disguised pointers to memory addresses in the heap (the only 
exceptions are so-called primitives, ``int``, ``double``, ``char``, etc). This means that
the arguments of factor graphs will always point to the same memory location in the heap.
This ensures that if two different factors share a variable at the beginning, they will always 
share it throughout sampling. Note that even if the *reference* is final, we can still modify
what is the *contents* of this memory location. Under the hood, this is what the MCMC sampler
does (for example for real variables, via ``variable.setValue()``.

#### Improving the factor via object and generic programming

We will now do a few things to improve our implementation.

**Static methods:** It is often convenient to access the logic inside a method (function associated
with an object) without always 
having to instantiate an object (new ClassName().function()). As a consequence, when 
possible it is good to pull the 
contents of method into a ``static`` function. See ``logDensity(double point, double rate)``
in ``Exponential.java``.
This way, it is possible to get the density using ``Exponential.logDensity(x,r)``, or even 
shorter, by adding ``import static bayonet.distributions.Exponential.*`` at the top of the 
file and then accessing the function with just ``logDensity(x,r)``.

**Flexible parameterization:** Distribution families are often parameterized using more 
than one type of parameterization.
For example, we would like the same Gamma class to support both shape-scale and shape-rate 
parameterizations. Here are the detailed steps of how this can be done:

1. Create an interface called ``Parameters``. In cases like here where an interface or class A is used mostly 
   in the context of another one B, it is a good idea to write A inside the file B.java as well.
   To do this, simply use ``public static interface B {}`` (``static`` here means that B does not
   have access to the fields in A; this is the best choice most of the time).
2. Create two classes ``ShapeScaleParameterization`` and ``ShapeRateParameterization`` 
   using the same strategy, one where you copy the shape and rate from earlier (keeping
   their annotations), and the second where you create two new fields with again the appropriate 
   annotations. Make all the fields final and create constructors. Also, make the two class implement Parameters.
3. Go back to the interface you have created, declare the signature of two methods, one that returns a rate,
   ``public double getRate();``
   one that return a shape.
4. Eclipse will now underline in red the two classes you have created. This is because they do not 
   respect the contract set by their interfaces. Clicking on the corresponding red icon
   in the margin, and selecting ``Add unimplemented methods`` will create stubs for the missing methods.
   Fill these methods: in one class, this will just return the identity, in the other, it does some simple
   transformations.
5. Finally, in the Gamma class, replace the shape and rate fields by an instance of Parameters, called ``parameters``. 
   Replace occurrences of ``rate.getValue()`` by ``parameters.getRate()``, and similarly for the shape. 
   Update the constructors.
6. Annotate the new field parameters with ``@FactorComponent``. This means that while parameters is not
   by itself a variable connected to this factor, it does contain variables connected to this factor 
   (those in the implementation annotated by ``@FactorArgument``. Note that this can be recursive (a 
   factor component can contain further factor components, etc).

**Generics:** The class as-is is full functional, but using it in method definitions is a bit cumbersome:
for example if we wanted to put an exponential prior on the rate of a gamma, we would need to write something 
like:

```java
@DefineFactor 
Gamma gamma = new Gamma(x);

@DefineFactor
Exponential exponential = Exponential.on(((RateParameterization) gamma.parameters).rate);
```

This is because in its current form, the compiler only knows that parameters is of type Parameter,
it does not known the precise type (``ShapeScaleParameterization`` or ``ShapeRateParameterization``).

To address this, we will use generics, which give us a way to declare which flavor of Gamma we use right when
we declare it. The end result will be being able to write:

```java
@DefineFactor 
Gamma<ShapeRateParameterization> gamma = new Gamma<ShapeRateParameterization>(x);

@DefineFactor
Exponential exponential = Exponential.on(gamma.parameters.rate);
```

To do this, proceed as follows:
1. Declare a generic type in the class declaration: ``public class Gamma<P extends Parameters>``. Here P is what is 
   called a *generic*. It is like a variable, but instead of holding a value, it holds a type (class or interface name).
   The part ``extends Parameters`` is called a type bound: it means that P is has to implement Parameters.
2. Change ``public final Parameters parameters;`` to ``public final P parameters;``.

**Method chaining:** This is a technique that again makes it just cleaner to  define models. This allows writing 
``Gamma.on(x).withShapeAndRate(y, z)``. See examples at the 
bottom of ``Exponential.java``.



#### Testing probabilistic programs

Writing probabilistic programs is tricky: subtle and not-so-subtle 
bugs can easily go unnoticed if one is not careful.

Blang offers some tools to help discover bugs in probabilistic programs.
The main one, ``CheckStationarity``, is based on the fact that forward sampling
is often  easier to implement than posterior sampling, and more  
importantly, is typically based on different ideas than posterior sampling.
Both types of sampling therefore share minimal
code with the posterior sampling code. So the hope is that they do not share
bugs either, so that discrepancy between the two methods can be detected.

Setting up this test is fairly straightforward:

```java
// Creates a model similar to the small example, but with no observations
// (forward simulation requires that there be no observations).
TestBlangSmallExample runner = new TestBlangSmallExample();

// Build an MCMC algorithm on this model (this was done implicitly in runner.run() earlier,
// but this time we will need a bit more control)
MCMCAlgorithm algo = runner.buildMCMCAlgorithm();

// Print a summary of the model to make sure the model is what we intend to test
// The syntax ${someString} in this output denotes a random variable named someString.
System.out.println("Summary of model");
System.out.println("----------------");
System.out.println(algo.model);

// The number of MCMC sweeps to perform at each iteration of the test (more info below)
algo.options.nMCMCSweeps = 10;

// Actual code for setting up the test itself
CheckStationarity check = new CheckStationarity();
check.setShowSampleSummaryStats(true);
System.out.println("Summary statistics of the samples");
System.out.println("---------------------------------");

// Here: 1000 is the number of test iterations (different than MCMC sweeps, see below)
//       0.05 is a p-value threshold
check.check(algo, 10000, 0.05);
```
<sub>From:[conifer.TestBlangSmallExample](src/test/java//conifer/TestBlangSmallExample.java)</sub>

Optional question: look at the information in ``CheckStationary.java``, and 
explain how it works and the theory behind it.

Note: this test can be quite expensive for large problems. Fortunately, most 
bugs can be detected on fairly small artificial datasets. So you should make your 
model/data you test against small enough so that you can frequently and quickly 
run test, but large enough so that you have *code coverage*, i.e. all parts of 
your code are needed to cover the examples being tested.

An important related principle 
is that you should strive to test against *all* small problem instances if possible.
When it is not possible, do as many as possible, focusing on corner and degenerate 
cases.

Also, make sure you test your test: commit your code, introduce one small bug, 
and see if your test can detect it.


To conclude this section, a second test, which makes sure that the randomness is
fixed, meaning that if you run two times the code with the same seed (as 
encoded in the ``Random`` object), then the exact same result should be obtained.

Very important for reproducibility and debugging.

```java
TestBlangSmallExample runner = new TestBlangSmallExample();
CheckFixedRandomness checkRand = new CheckFixedRandomness(runner.factory);
checkRand.check(runner);
```
<sub>From:[conifer.TestBlangSmallExample](src/test/java//conifer/TestBlangSmallExample.java)</sub>

Bayesian Phylogenetic inference
-------------------------------

We will now apply this framework to a more complex problem, Bayesian phylogenetic
inference. We start with the core of the likelihood calculation, the matrix exponential.

#### Implementing the matrix exponential



Fill in ``marginalTransitionProbability()``.

Use the diagonalization method covered in class, using 
the eigen-decomposition functionalities provided by EJML.

<sub>From:[conifer.ctmc.CTMC](src/main/java//conifer/ctmc/CTMC.java)</sub>

Fill in ``stationaryDistribution()``.

Show how you can reduce the problem of finding the stationary distribution
of a CTMC to a eigenvalue problem. Again, solve this problem using
EJML.

<sub>From:[conifer.ctmc.CTMC](src/main/java//conifer/ctmc/CTMC.java)</sub>

#### Computing the likelihood on a tree




<sub>From:[conifer.factors.TreeLikelihood](src/main/java//conifer/factors/TreeLikelihood.java)</sub>

