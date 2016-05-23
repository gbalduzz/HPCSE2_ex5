Giovanni Balduzzi's exercise 5 for HPCSE2. Work in progress.

Branch description:
Naive: compute expansion from scratches for every node
Naive_e2e: it's only purpose is to test the e2e results with against a "stupid" e2e kernel that uses std::complex
e2e: the interesting branch: I rewrote the sum inside the expansion translation formula as a sum over the powers of z.
     Therefore each power  z**p is computed once and then it contributes to all the coefficients of order l>=p.
     Unfortunately I used some ugly template recursion for unrolling the loop. I wonder if learning m4 would provide better time to code and readability.
buffered_e2p: I tried to implement the buffering as shown in the slides. Unfortunately the code does note vectorize with gcc + openmp. 
