# Frequently Asked Questions

1. Why am I getting segmentation faults after edge flipping?

Edge flipping changes the topology of the underlying mesh. The previously 
computed orientation is now invalid. You must recompute the orientation of the 
simplicial complex using the following function:
```
compute_orientation(mesh);
```
In the future we hope to have on-the-fly computation of orientation. But we are
not quite there yet!