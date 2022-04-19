# Skeletal Representations

Skeletal representations (s-reps) have proved to be a robust model for shape analysis. This extension allows for the creation, refinement, and viewing of these s-reps within [3D Slicer](https://www.slicer.org/). It has also been incorporated into [SlicerSALT](https://salt.slicer.org/).

This extension introduces a new MRML data node interface, `vtkMRMLSRepNode`, and an implementation, `vtkMRMLEllipticalSRepNode`. Additionally, other MRML nodes are added for display and storage purposes. `vtkMRMLSRepNode` allow for programmatic editing of s-reps. These new MRML nodes can be saved and loaded to `.srep.json` files via Slicer's normal saving and loading mechanisms.

## Modules

| Name | Description |
|------|-------------|
| [SRep](SRep) | The base module for s-reps. This module introduces the MRML nodes, and has generic facilities for display. |
| [SRepCreator](SRepCreator) | Allows creating initial fit s-reps via flowing a model to an ellipse, generating an s-rep for that ellipse, and then backflowing the s-rep via thin plate spline to fit the original model. Essentially the SRepInitializer from the previous implementation. |
| [SRepRefinement](SRepRefinement) | Allows optimizing an s-rep fit via min_newuoa optimization. Essentially the SRepRefiner from the previous implementation. |

## Resources

To learn more about Slicer, SlicerSALT, and Slicer extensions, check out the following resources.

 - https://slicer.readthedocs.io/en/latest/
 - https://salt.slicer.org/
 - https://slicer.readthedocs.io/en/latest/user_guide/extensions_manager.html

To learn more about s-reps, see the following papers.

 - Liu, Z., Hong, J., Vicory, J., Damon, J. N., & Pizer, S. M. (2021). Fitting unbranching skeletal structures to objects. _Medical Image Analysis_, 70, 102020.

## History

There was a previous implementation of this extension that pioneered the s-rep work inside of SlicerSALT. That implementation did not have any unified classes for working with SReps. This implementation sought to bring a common way of reading, writing, and interacting with s-reps. With the creation of the new MRML nodes to do this, the other modules needed significant updates to works with them. Hence the new implementation.

During development, the two s-rep implementations lived in parallel, so the new could be tested against the old. This is the reason that none of the modules kept the same name, as 3D Slicer does not allow modules with identical names.

## License

This software is licensed under the terms of the [Apache Licence Version 2.0](LICENSE).

