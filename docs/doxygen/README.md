
## Files used for the documentation

The following files are used for creating the documentation:

 - **README.md**: This file.
 - **config.doxygen**: This file is setup with the options we want to configure in doxygen.
 - **footerfile.html, headerfile.html**: These files are used as the header and footer of your html documentation. I have specified the icon in the header here.
 - **forest_icon.png**: The icon for the aplication. Not necessary at all.
 - **MyLayout.xml**: I specify here if some special layout should be used.
 - **mainpage.dox, diffusion.dox, input.dox, multiphysics.dox, neutronics.dox, output.dox, results.dox, thermalhydraulics.dox, transport.dox**: Different files to explain the problem before showing the code.
 - **science.bib**: The bibliography file to cite proper references in the explanations before.
 - **tutorial.dox, test1.dox test2.dox**: Here we show how to solve the problem without going too much inside the code.


## How to build the documentation

In order to construct the documentation with Doxygen we type

```console
doxygen config.doxygen
```

We can then visualize the generated html, for instance with firefox, by typing

```console
firefox html/index.html &
```


