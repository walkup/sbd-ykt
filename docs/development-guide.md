# Development Guide

## For source code writing
### Macro protection format
`SBD_<FILE_PATH>`. For example: `SBD_SBD_H` for `sbd/sbd.h` file and `SBD_CHEMISTRY_FOO_BAR_H` for `sbd/chemistry/foo_bar.h` file. A complete example looks like

```cpp
#ifndef SBD_CHEMISTRY_FOO_BAR_H
#define SBD_CHEMISTRY_FOO_BAR_H


// code from here


#endif /* ifndef SBD_CHEMISTRY_FOO_BAR_H */
```

### Class structure format
```cpp
class FooBar {
public:
  PublicFunc1();
  PublicFunc2();

  Type public_var1;
  Type public_var2;

private:
  Type private_var1_;
  Type private_var2_;


  PrivateFun1_();
  PrivateFun2_();
};
```

### Indention
2 spaces.

### Maximal line length
Try to keep it shorter than 80. Please no longer than 90 (not including comments)!

### Internal file include
Always include file from root path, for example, in the file `sbd/foo/bar.h`
```cpp
#include "sbd/sbd.h"
```

### MPI related code writing
*Under updating ...*

## For code documenting
### Write the document
We use [Doxygen](https://www.doxygen.org/index.html) to document the code. Documenting is boring but it is extremely important for future reference. So please treat the document as careful as the source code. Here is an example about the document of a file named `foo.h`.

```
//// This file is a part of GraceQ/tensor
/**
@file foo.h
@brief A brief description of this file.
*/
#ifndef SBD_FOO_H
#define SBD_FOO_H


namespace gqten {


/**
Short description of class Foo. Detailed description of class Foo.

@tparam T Description of T.
*/
template <typename T>
class Foo {
public:
  /**
  Short description of member function PublicFunc1. Detailed description if it is needed.

  @param[in, out] var1 Description of var1.

  @note Write down some notes for this function.

  #### Examples
  @code
  Foo foo();
  foo.PublicFunc1(bar);
  // Possible behavior
  @ endcode   # remove the space between "@" and "endcode" in the real source code
  */
  PublicFunc1(Type var1);

  /// For small function, one-line documentation is ok.
  PublicFunc2();

  /// One-line documentation of public_var1.
  Type public_var1;

  /// One-line documentation of public_var2.
  Type public_var2;

private:

  /**
  Private members also need document. Please treat them as well as public members.
  */
  Type private_var1_;

  /// Document of private_var2_.
  Type private_var2_;

  /// Document of PrivateFun1_.
  PrivateFun1_();

  /// Document of PrivateFun2_.
  PrivateFun2_();
};


}
#endif /* ifndef SBD_FOO_H */
```

### Build the document
```sh
> make doc/build    # generate document to build_doxy/ folder

> make doc/watch    # Watch the modification of the files and rebuild the document at real time (need fswatch installed)
```


## For testing
### Write the test
Writing test is also very boring but again extremely important. Tests in the lower level library let you keep away from subtle bugs which you definitely do not want to meet when you build your applications. So please treat them as careful as the source code, if you do not want your excited *new physics* actually is a *new bug*.

### Build and run the tests
See the ./sample/sample_based_diagonalization/README.md

## For Git usage
### commit message
The commit message should the clear information about what is included in this commit. Please avoid something like `"A small fix"`(what is fixed?) or `"Performance optimization"`(for which function?).
