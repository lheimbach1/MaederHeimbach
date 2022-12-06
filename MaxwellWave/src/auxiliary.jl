@doc raw"""
    save_array(Aname,A)

Saves the array A in a binary file with name Aname.
```
fname = string(Aname,".bin")
out = open(fname,"w"); write(out,A); close(out)
```
"""
function save_array(Aname,A)
    fname = string(Aname,".bin")
    out = open(fname,"w"); write(out,A); close(out)
end

@doc raw"""
    load_array(Aname,A)

Loades the array A from a binary file with name Aname.
```
fname = string(Aname,".bin")
fid=open(fname,"r"); read!(fid,A); close(fid)
```
"""
function load_array(Aname,A)
    fname = string(Aname,".bin")
    fid=open(fname,"r"); read!(fid,A); close(fid)
end
