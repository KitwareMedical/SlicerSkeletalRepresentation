Nick and Zhiyuan's sreps have different ordering of the srep points despite having the same structure

Zhiyuan created a program to reorder the points in Nick's srep so that Zhiyuan's refinement could work

Nick modified that program so that it kept track of the reordering changes that occurred and stored a list called mapping.pkl

Nick created a program that takes mapping.pkl and a re-ordered srep to produce an srep that is back to nick's ordering

Now we can do:
1. Re-order nick's srep
2. Refine the re-ordered srep
3. Convert the refined srep back to nick's original ordering