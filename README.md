# GSpyrou_package

## **Viterbi and ForwardBackward algorithms for Symmetric Continuous-time Markov-Chain HMM**

This R package contains two main algorithms for generic Continuous-time Markov Chain hidden Markov Models and its written in C++

**1) Implementation of ForwardBackward Algorithm on a Symmetric Continuous-time Markov Chain HMM:**

 Function that performs the ForwardBackward algorithm,given likelihood of observed data,timestamp of each observation and a parameter  theta.

**2) Implementation of Viterbi Algorithm on a Symmetric Continuous-time Markov Chain HMM using Rcpp:**

 Function that implements a Viterbi algorithm, given likelihood of observed data,timestamp of each observation and a parameter theta.


You can download the package in your computer or just import it directly into your R environment by writing the following command:
devtools::install_github("gpsyrou/GSpyrou_package")

**References**

https://en.wikipedia.org/wiki/Viterbi_algorithm

https://en.wikipedia.org/wiki/Forward-backward_algorithm
