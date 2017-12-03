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



2. I'm getting errors like 
```
error: expected expression
        size_t nverts = mesh.size<1>();
        								 ^
```
, what gives?


Sometimes metatemplate programming complicates things. In this case, the
compiler is having a hard time understanding what `.size<1>()` means. We can 
disambiguate by telling it that it refers to a function template call.
```
size_t nverts = mesh.template size<1>();
```
Read more about it here [ http://en.cppreference.com/w/cpp/language/dependent_name ]