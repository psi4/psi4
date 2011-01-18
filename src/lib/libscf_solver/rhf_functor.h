#ifndef __rhf_functor_h__
#define __rhf_functor_h__

namespace psi { namespace scf {

class SOERI {
    const shared_ptr<Matrix> &D_;
    shared_ptr<Matrix> &G_;

    inline int integral_type(int i, int j, int k, int l)
    {
        int type;

        if (i == j && i == k && i == l)     // (ij|kl)  (11|11)
            type = 1;
        else if (i == j && k == l && i > k) // (ij|kl)  (22|11)
            type = 2;
        else if (i == j && i == k && i > l) // (ij|kl)  (22|21)
            type = 3;
        else if (j == k && j == l && i > j) // (ij|kl)  (21|11)
            type = 4;
        else if (i == k && j == l)          // (ij|kl)  (21|21)
            type = 5;
        else if (i == j)                    // (ij|kl)  (33|21)
            type = 6;
        else if (j >  k && k == l)          // (ij|kl)  (32|11)
            type = 7;
        else if (k == l)                    // (ij|kl)  (31|22)
            type = 8;
        else if (i == k)                    // (ij|kl)  (32|31)
            type = 9;
        else if (j == k)                    // (ij|kl)  (32|21)
            type = 10;
        else if (j == l)                    // (ij|kl)  (31|21)
            type = 11;
        else if (j >  k)                    // (ij|kl)  (43|21)
            type = 12;
        else if (j >  l)                    // (ij|kl)  (42|31)
            type = 13;
        else                                // (ij|kl)  (41|32)
            type = 14;

        return type;
    }

public:
    SOERI(const shared_ptr<Matrix>& D, shared_ptr<Matrix>& G)
        : D_(D), G_(G), counter(0)
    { }

    size_t counter;
    void operator()(int i, int j, int k, int l, int is, int ii, int js, int jj, int ks, int kk, int ls, int ll, double value) {
    double temp1=0.0, temp2=0.0, temp3=0.0, temp4=0.0, temp5=0.0, temp6=0.0;
    int itype;

        fprintf(outfile, " (%2d %2d %2d %2d) = %14.8lf\n", i, j, k, l, value);

            itype = integral_type(i, j, k, l);
            switch(itype) {
                case 1:
                    temp1 = D_->get(is, ii, ii) * value;

                    G_->add(is, ii, ii, temp1);
                    break;

                    case 2:
                    temp1 = D_->get(ks, kk, kk) * 2.0 * value;
                    temp2 = 0.0;
                    temp3 = D_->get(is, ii, ii) * 2.0 * value;

                    G_->add(is, ii, ii, temp1);
                    G_->add(ks, kk, kk, temp3);

                    if (is == ks) {
                        temp2 = D_->get(is, ii, kk) * value;
                        G_->add(is, ii, kk, -temp2);
                        G_->add(is, kk, ii, -temp2);
                    }
                    break;

                    case 3:
                    temp1 = temp2 = 0.0;
                    if (is == ls) {
                        temp1 = D_->get(is, ii, ii) * value;
                        temp2 = D_->get(is, ii, ll) * value * 2.0;

                        G_->add(is, ii, ll, temp1);
                        G_->add(is, ll, ii, temp1);
                        G_->add(is, ii, ii, temp2);
                    }
                    break;

                    case 4:
                    temp1 = temp2 = 0.0;
                    if (is == js) {
                        temp1 = D_->get(js, jj, jj) * value;
                        temp2 = D_->get(is, ii, jj) * value * 2.0;

                        G_->add(is, ii, jj, temp1);
                        G_->add(is, jj, ii, temp1);
                        G_->add(js, jj, jj, temp2);
                    }
                    break;

                    case 5:
                    if (is == js) {
                        temp1 = D_->get(is, ii, jj) * value * 3.0;
                        G_->add(is, ii, jj, temp1);
                        G_->add(is, jj, ii, temp1);
                    }

                    temp2 = D_->get(is, ii, ii) * value;
                    temp3 = D_->get(js, jj, jj) * value;
                    G_->add(js, jj, jj, -temp2);
                    G_->add(is, ii, ii, -temp3);
                    break;

                    case 6:
                    temp1 = temp2 = temp3 = temp4 = 0.0;
                    if (ks == ls)
                        temp1 = D_->get(ks, kk, ll) * value * 4.0;
                    if (is == ls)
                        temp2 = D_->get(is, ii, ll) * value;
                    temp3 = D_->get(is, ii, ii) * value * 2.0;
                    if (is == ks)
                        temp4 = D_->get(is, ii, kk) * value;

                    G_->add(is, ii, ii, temp1);
                    if (is == ks) {
                        G_->add(is, ii, kk, -temp2);
                        G_->add(is, kk, ii, -temp2);
                    }
                    if (ks == ls) {
                        G_->add(ks, kk, ll, temp3);
                        G_->add(ks, ll, kk, temp3);
                    }
                    if (is == ls) {
                        G_->add(is, ii, ll, -temp4);
                        G_->add(is, ll, ii, -temp4);
                    }
                    break;

                    case 7:
                    temp1 = temp2 = temp3 = temp4 = 0.0;
                    if (is == js)
                        temp1 = D_->get(is, ii, jj) * value * 4.0;
                    if (js == ks)
                        temp2 = D_->get(js, jj, kk) * value;
                    if (is == ks)
                        temp3 = D_->get(is, ii, kk) * value;
                    temp4 = D_->get(ks, kk, kk) * value * 2.0;

                    G_->add(ks, kk, kk, temp1);
                    if (is == ks) {
                        G_->add(is, ii, kk, -temp2);
                        G_->add(is, kk, ii, -temp2);
                    }
                    if (js == ks) {
                        G_->add(js, jj, kk, -temp3);
                        G_->add(js, kk, jj, -temp3);
                    }
                    if (is == js) {
                        G_->add(is, ii, jj, temp4);
                        G_->add(is, jj, ii, temp4);
                    }
                    break;

                    case 8:
                    temp1 = temp2 = temp3 = temp4 = 0.0;
                    temp1 = D_->get(ks, kk, kk) * value * 2.0;
                    if (is == js)
                        temp2 = D_->get(is, ii, jj) * value * 4.0;
                    if (js == ks)
                        temp3 = D_->get(js, jj, kk) * value;
                    if (is == ks)
                        temp4 = D_->get(is, ii, kk) * value;

                    if (is == js) {
                        G_->add(is, ii, jj, temp1);
                        G_->add(is, jj, ii, temp1);
                    }
                    G_->add(ks, kk, kk, temp2);
                    if (is == ks) {
                        G_->add(is, ii, kk, -temp3);
                        G_->add(is, kk, ii, -temp3);
                    }
                    if (js == ks) {
                        G_->add(js, jj, kk, -temp4);
                        G_->add(js, kk, jj, -temp4);
                    }
                    break;

                    case 9:
                    temp1 = temp2 = temp3 = temp4 = 0.0;
                    if (is == ls)
                        temp1 = D_->get(is, ii, ll) * value * 3.0;
                    if (is == js)
                        temp2 = D_->get(is, ii, jj) * value * 3.0;
                    if (js == ls)
                        temp3 = D_->get(js, jj, ll) * value * 2.0;
                    temp4 = D_->get(is, ii, ii) * value;

                    if (is == js) {
                        G_->add(is, ii, jj, temp1);
                        G_->add(is, jj, ii, temp1);
                    }
                    if (is == ls) {
                        G_->add(is, ii, ll, temp2);
                        G_->add(is, ll, ii, temp2);
                    }
                    G_->add(is, ii, ii, -temp3);
                    if (js == ls) {
                        G_->add(js, jj, ll, -temp4);
                        G_->add(js, ll, jj, -temp4);
                    }
                    break;

                    case 10:
                    temp1 = temp2 = temp3 = temp4 = 0.0;
                    if (js == ls)
                        temp1 = D_->get(js, jj, ll) * value * 3.0;
                    if (is == js)
                        temp2 = D_->get(is, ii, jj) * value * 3.0;
                    temp3 = D_->get(js, jj, jj) * value;
                    if (is == ls)
                        temp4 = D_->get(is, ii, ll) * value * 2.0;

                    if (is == js) {
                        G_->add(is, ii, jj, temp1);
                        G_->add(is, jj, ii, temp1);
                    }
                    if (js == ls) {
                        G_->add(js, jj, ll, temp2);
                        G_->add(js, ll, jj, temp2);
                    }
                    if (is == ls) {
                        G_->add(is, ii, ll, -temp3);
                        G_->add(is, ll, ii, -temp3);
                    }
                    G_->add(js, jj, jj, -temp4);
                    break;

                    case 11:
                    temp1 = temp2 = temp3 = temp4 = 0.0;
                    if (ks == js)
                        temp1 = D_->get(ks, kk, jj) * value * 3.0;
                    if (is == js)
                        temp2 = D_->get(is, ii, jj) * value * 3.0;
                    temp3 = D_->get(js, jj, jj) * value;
                    if (is == ks)
                        temp4 = D_->get(is, ii, kk) * value * 2.0;

                    if (is == js) {
                        G_->add(is, ii, jj, temp1);
                        G_->add(is, jj, ii, temp1);
                    }
                    if (ks == js) {
                        G_->add(ks, kk, jj, temp2);
                        G_->add(ks, jj, kk, temp2);
                    }
                    if (is == ks) {
                        G_->add(is, ii, kk, -temp3);
                        G_->add(is, kk, ii, -temp3);
                    }
                    G_->add(js, jj, jj, -temp4);
                    break;

                    case 12:
                    case 13:
                    case 14:
                    temp1 = temp2 = temp3 = temp4 = temp5 = temp6 = 0.0;
                    if (ks == ls)
                        temp1 = D_->get(ks, kk, ll) * value * 4.0;
                    if (is == js)
                        temp2 = D_->get(is, ii, jj) * value * 4.0;
                    if (js == ls)
                        temp3 = D_->get(js, jj, ll) * value;
                    if (is == ks)
                        temp4 = D_->get(is, ii, kk) * value;
                    if (js == ks)
                        temp5 = D_->get(js, jj, kk) * value;
                    if (is == ls)
                        temp6 = D_->get(is, ii, ll) * value;

                    if (is == js) {
                        G_->add(is, ii, jj, temp1);
                        G_->add(is, jj, ii, temp1);
                    }
                    if (ks == ls) {
                        G_->add(ks, kk, ll, temp2);
                        G_->add(ks, ll, kk, temp2);
                    }
                    if (is == ks) {
                        G_->add(is, ii, kk, -temp3);
                        G_->add(is, kk, ii, -temp3);
                    }
                    if (js == ls) {
                        G_->add(js, jj, ll, -temp4);
                        G_->add(js, ll, jj, -temp4);
                    }
                    if (is == ls) {
                        G_->add(is, ii, ll, -temp5);
                        G_->add(is, ll, ii, -temp5);
                    }
                    if (js == ks) {
                        G_->add(js, jj, kk, -temp6);
                        G_->add(js, kk, jj, -temp6);
                    }
                    break;
            };
            counter++;
        }
};

} }

#endif

