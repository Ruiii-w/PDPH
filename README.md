# PDPH
PDPH: A Performant Dynamic Perfect Hash Index Using Hybrid
Strategies for Storage Organizations

We introduce PDPH, a performant dynamic perfect
hash index with hybrid storage organization that hugely reduces
rehashing overhead. (1) PDPH applies loose bucket organizations
with a strict bucket isolation policy, which are design decisions that
other methods have overlooked within the design space. The loose
bucket organization eliminates redundant rehashing attempts, and
the strict bucket isolation policy prevents access to unused records,
which mitigates the above amplification issue during rehashing
processes at the cost of temporary space overhead. (2) PDPH pro-
poses multiple-level hybrid structures for buckets with different
utilization, respectively. It exploits the unevenness of the distribu-
tion between hash buckets with different utilization to handle the
potential memory inefficiency caused by loose bucket organizations
while preserving the query performance.
