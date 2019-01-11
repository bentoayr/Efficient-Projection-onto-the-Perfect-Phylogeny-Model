# Efficient Projection onto the Perfect Phylogeny Model


This code implements a generalization of the algorithm described in the paper "[Efficient Projection onto the Perfect Phylogeny Model](https://arxiv.org/pdf/1811.01129.pdf)".

Please cite this code as:

@inproceedings{projection-onto-PPM,
  title={Efficient Projection onto the Perfect Phylogeny Model},
  author={Jia, Bei and Ray, Surjyendu  and Safavi, Sam and Bento, Jose},
  booktitle={Advances in Neural Information Processing Systems},
  year={2018}
}

To compile the code just do gcc project_onto_PPM.c 

To run the code just do ./a.out InputData OutputData 0     or   ./a.out InputData OutputData 1



The dependencies of the different functions in the code is described by the following diagram.

![alt text](https://raw.githubusercontent.com/bentoayr/Efficient-Projection-onto-the-Perfect-Phylogeny-Model/master/pic/cflow0.png)
