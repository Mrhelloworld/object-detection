package com.atex.csce482.objectdetection;

/**
 * Created by Chuong on 3/21/2015.
 */
public class Kernel {
    static double[][] kernel;
    static int kval;

    public Kernel(int k) {

        kval = k;
        kernel = new double[kval][kval];
        return;
    }
}
