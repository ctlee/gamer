
// %typemap(out) (std::string) { $result = PyString_FromString($1.c_str());}
