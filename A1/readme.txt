Approach 1: Computing all elements then Outputting Upper Triangle
Approach 2: Computing only Upper Triangle

1:Idea
    Divide the output matrix into threads row wise, compute all elements and only output upper triangle
1:Why Attempted
    Obvious as to how work is to be divided among threads
1:Results
    Scaled decently
    for O(2^10) non zero input elements, time taken was 0.12s
    for O(2^24) non zero input elements, time taken was 38m04.45s
    Better than O(n^3) scaling expected for matrix multiplacation
1:Drawback
    Features a large amount of redundant computation, i. e. the useful work done by threads is highly skewed.

2:Idea
    Divide the output matrix into threads in such a way that each thread does an equal amount of work
2:Why Attempted
    Each thread has an equal amount of useful work so maximum parallelism is exploited
2:Results
    O(n^1.5) type scaling - far better than O(n^3) expected by naive implementation
    for O(2^10) non zero input elements, time taken was 0.01s
    for O(2^24) non zero input elements, time taken was 8.02s
2:Drawback
    Tasks are of different sizes so maximum parallelism is not exploited


Final Scalability Analysis

Non-0 Input Blocks, Non-0 Output Blocks, Runtime (s)
93(Approx 2^7),210,0.013
4095(Approx 2^12),4095,1.095
19923(Approx 2^14),200024,8.023

Note: Testing done on local machine as css was not working consistently.