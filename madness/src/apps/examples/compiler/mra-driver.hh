#ifndef MRA_DRIVER_HH
# define MRA_DRIVER_HH
# include <string>
# include <map>
# include <vector>
# include <sstream>
# include "mra-parser.hh"

// Tell Flex the lexer's prototype ...
# define YY_DECL                                      \
  yy::mra_parser::token_type                       \
  yylex (yy::mra_parser::semantic_type* yylval,      \
         yy::mra_parser::location_type* yylloc,      \
         mra_driver& driver)
// ... and declare it for the parser's sake.
YY_DECL;

class Exp;

extern std::map<std::string, bool> dectab;        // keeps track of CXX declared symbols

std::ostream& operator<<(std::ostream& s, const Exp& e);

class Exp {
public:
    enum Type {
        INTEGER, 
        REAL, 
        COMPLEX, 
        REALTEN, 
        COMPLEXTEN, 
        REALFUN, 
        COMPLEXFUN, 
        STRING, 
        COORD,
        UNKNOWN
    };

    static const char* type_name(Type t) {
        static const char* name[] = {"Z ", "R ", "C ", "RT", "CT", "RF", "CF", "S ", "Rd", "U "};
        return name[t];
    }

    static const char* cxx_type_name(Type t) {
        static const char* name[] = {"long", "double", "double_complex", "real_tensor", "complex_tensor", "real_function", "complex_function", "std::string", "coord", "unknown"};
        return name[t];
    }
    
    typedef std::vector<Exp*>::iterator iterator;
    typedef std::vector<Exp*>::const_iterator const_iterator;

private:
    std::string s;       // Token or operator
    Type t;              // Result type
    std::vector<Exp*> c; // Children in order

public:

    // Leaf node
    Exp(const std::string& s, Type t=UNKNOWN) 
        : s(s), t(t), c()
    {}

    // Unary op
    Exp(const std::string& s, Exp* left, Type t=UNKNOWN) 
        : s(s), t(t), c()
    {
        c.push_back(left);
        if (is_unknown(t)) this->t = unaryoptype(s,left->type());
    }

    // Binary op
    Exp(const std::string& s, Exp* left, Exp*right, Type t=UNKNOWN)
        : s(s), t(t), c()
    {
        c.push_back(left);
        c.push_back(right);
        if (is_unknown(t)) this->t = binaryoptype(s,left->type(),right->type());
    }

    // Ternary op
    Exp(const std::string& s, Exp* left, Exp* middle, Exp* right, Type t=UNKNOWN)
        : s(s), t(t), c()
    {
        c.push_back(left);
        c.push_back(middle);
        c.push_back(right);
    }

    // Quaternary op
    Exp(const std::string& s, Exp* left, Exp* middle, Exp* right, Exp* far, Type t=UNKNOWN)
        : s(s), t(t), c()
    {
        c.push_back(left);
        c.push_back(middle);
        c.push_back(right);
        c.push_back(far);
    }

    Exp* add(Exp* child) 
    {
        c.push_back(child);
        return child;
    }

    Exp* add(const std::string& str, Type t=UNKNOWN) 
    {
        return add(new Exp(str,t));
    }

    void print_tree(std::ostream& file, int depth=0) const {
        for (int i=0; i<depth; i++) file << "        ";
        file << *this << std::endl;
        for (const_iterator iter = c.begin(); 
             iter != c.end(); 
             ++iter) (*iter)->print_tree(file,depth+1);
        
    }

    bool match_sub_tree(const Exp* t) const {
        // t points to a template tree that we want to
        // match against this sub-tree. either exact match
        // or wildcard indicated by null string in t
        
        if (t->str() == "" || (t->str() == str() && t->c.size() <= c.size())) {
            const_iterator it = t->begin();
            const_iterator ti =    begin();
            while (it != t->end()) {
                if (! (*ti)->match_sub_tree(*it)) return false;
                ++it;
                ++ti;
            }
            return true;
        }
        return false;
    }

    void regenerate(std::ostream& file) const {
        // varlist, arglist, paramlist can be 1 or more entries
        if (s == "arglist" || s=="varlist" || s=="paramlist") {
            const_iterator iter = c.begin(); 
            while (iter != c.end()) {
                (*iter)->regenerate(file);
                ++iter;
                if (iter == c.end()) {
                    break;
                }
                else {
                    file << ", ";
                }
            }
        }
        else if (c.size() == 0) {
            file << s;
        }
        else if (c.size() == 1) {
            // ()       (exp)
            // -        - exp
            // not      not exp
            // greek    greek name
            // print    print arglist
            // plot     plot  exp
            // F        F variable
            // while    while exp
            // if       if exp
            // elif     elif exp
            // trace    trace exp

            if      (s == "()") file << "(";
            else if (s == "-") file << "-";
            else if (s == "greek");
            else if (s == "trace") file << "<"; 
            else file << s << " ";

            c[0]->regenerate(file);

            if (s == "()") file << ")";
            else if (s == "trace") file << ">"; 
        }
        else if (c.size() == 2) {
            // decl     varlist in field
            // call     variable(arglist)
            // mathfun  funcname(exp)
            // braket   <a|b>
            // =, +, -, *, _, ^    exp op exp
            // del

            if (s == "decl") {
                c[0]->regenerate(file);
                file << "in ";
                c[1]->regenerate(file);
            }
            else if (s == "call") {
                c[0]->regenerate(file);
                file << "(";
                c[1]->regenerate(file);
                file << ")";
            }
            else if (s == "mathfun") {
                if      (c[0]->str() == "abs") file << "|";
                else if (c[0]->str() == "norm") file << "||";
                else if (c[0]->str() == "normf") file << "||";
                else {
                    c[0]->regenerate(file);
                    file << "(";
                }
                c[1]->regenerate(file);
                if      (c[0]->str() == "abs") file << "|";
                else if (c[0]->str() == "norm") file << "||";
                else if (c[0]->str() == "normf") file << "||_F";
                else file << ")";
            }
            else if (s == "braket") {
                file << "<";
                c[0]->regenerate(file);
                file << "|";
                c[1]->regenerate(file);
                file << ">";
            }
            else if (s == "del") {
                file << "del_";
                c[0]->regenerate(file);
                file << " ";
                c[1]->regenerate(file);
            }
            else if (s=="_" || s=="^") {
                c[0]->regenerate(file);
                file << s;
                c[1]->regenerate(file);
            }
            else if (s == "lt") {
                c[0]->regenerate(file);
                file << " < ";
                c[1]->regenerate(file);
            }
            else if (s == "gt") {
                c[0]->regenerate(file);
                file << " < ";
                c[1]->regenerate(file);
            }
            else if (s == "le") {
                c[0]->regenerate(file);
                file << " <= ";
                c[1]->regenerate(file);
            }
            else if (s == "ge") {
                c[0]->regenerate(file);
                file << " >= ";
                c[1]->regenerate(file);
            }
            else if (s == "eq") {
                c[0]->regenerate(file);
                file << " == ";
                c[1]->regenerate(file);
            }
            else if (s == "ne") {
                c[0]->regenerate(file);
                file << " != ";
                c[1]->regenerate(file);
            }
            else {
                c[0]->regenerate(file);
                file << " " << s << " ";
                c[1]->regenerate(file);
            }
        }
        else if (c.size() == 3) {
            // def      variable = paramlist -> exp
            // domain   lo, hi, D
            // ^_       x^y_z

            if (s == "def") {
                c[0]->regenerate(file);
                file << " = ";
                c[1]->regenerate(file);
                file << " -> ";
                c[2]->regenerate(file);
            }
            else if (s == "domain") {
                file << "Omega = [";
                c[0]->regenerate(file);
                file << ", ";
                c[1]->regenerate(file);
                file << "]^";
                c[2]->regenerate(file);
            }
            else if (s == "^_") {
                c[0]->regenerate(file);
                file << "^";
                c[1]->regenerate(file);
                file << "_";
                c[2]->regenerate(file);
            }
            else if (s == "for") {
                file << "for ";
                c[0]->regenerate(file);
                file << " in [";
                c[1]->regenerate(file);
                file << ",";
                c[2]->regenerate(file);
                file << "]";
            }
            else {
                std::cout << "Trying to generate unknown expression with 3 arguments?\n";
                this->print_tree(std::cout);
            }
        }
        else if (s == "sum") {
            file << "sum_";
            c[0]->regenerate(file);
            file << "=";
            c[1]->regenerate(file);
            file << "^";
            c[2]->regenerate(file);
            file << " ";
            c[3]->regenerate(file);
        }
        else {
            std::cout << "Trying to generate expression with more than 4 arguments?\n";
            this->print_tree(std::cout);
        }
    }

    void generate_tex(std::ostream& file) const {
        // varlist, arglist, paramlist can be 1 or more entries
        if (s == "arglist" || s=="varlist" || s=="paramlist") {
            const_iterator iter = c.begin(); 
            while (iter != c.end()) {
                (*iter)->generate_tex(file);
                ++iter;
                if (iter == c.end()) {
                    break;
                }
                else {
                    file << ", ";
                }
            }
        }
        else if (c.size() == 0) {
            if (s=="Let" || s=="In" || s=="End" || s=="else" || s=="end" || s=="true" || s=="false" || s=="break") {
                file << "\\textsf{" << s <<"}\\ ";
            }
            else if (is_string(t)) {
                file << "\\verb+" << s << "+";
            }
            else {
                file << s;
            }
        }
        else if (c.size() == 1) {
            // ()       (exp)
            // {}       {exp}
            // -        - exp
            // not      not exp
            // F        F variable
            // greek    greek name
            // print    print arglist
            // if       if exp
            // elif     elif exp
            // while    while exp
            // plot     plot  exp
            // lap
            // trace    trace exp
            
            if      (s == "()") file << "\\left(";
            else if (s == "{}") file << "{";
            else if (s == "-") file << "-";
            else if (s == "F") file << "\\mathcal{F}\\ ";
            else if (s == "greek") file << "\\";
            else if (s == "lap") file << " \\nabla^2 ";
            else if (s == "trace") file << " \\langle ";
            else file << "\\textsf{" << s << "}\\ ";
            
            c[0]->generate_tex(file);
            
            if      (s == "{}") file << "}";
            else if (s == "()") file << "\\right)";
            else if (s == "trace") file << " \\rangle ";
        }
        else if (c.size() == 2) {
            // decl     varlist in field
            // call     variable(arglist)
            // mathfun  funcname(exp)
            // braket   <a|b>
            // =, +, -, *, _, ^    exp op exp
            // del
            // bsh

            if (s == "decl") {
                c[0]->generate_tex(file);
                file << " \\in ";
                c[1]->generate_tex(file);
            }
            else if (s == "call") {
                c[0]->generate_tex(file);
                file << "\\left(";
                c[1]->generate_tex(file);
                file << "\\right)";
            }
            else if (s == "mathfun") {
                file << " \\";
                c[0]->generate_tex(file);
                const std::string cs = c[0]->str();
                if (cs=="sqrt" || cs=="abs" || cs=="norm" || cs=="normf") file << "{";
                else file << "\\left(";
                c[1]->generate_tex(file);
                if (cs=="sqrt" || cs=="abs" || cs=="norm" || cs=="normf") file << "}";
                else file << "\\right)";
            }
            else if (s == "braket") {
                file << " \\langle ";
                c[0]->generate_tex(file);
                file << " | ";
                c[1]->generate_tex(file);
                file << " \\rangle ";
            }
            else if (s == "bsh") {
                if (c[0]->str() == "0.0") {
                    file << " \\nabla^{-2} \\left(";
                    c[1]->generate_tex(file);
                    file << "\\right) ";
                }
                else {
                    file << " \\left(";
                    c[0]->generate_tex(file);
                    file << " - \\nabla^2 \\right)^{-1} \\left(";
                    c[1]->generate_tex(file);
                    file << "\\right) ";
                }
            }
            else if (s == "=") {
                c[0]->generate_tex(file);
                file << " \\ = \\ ";
                c[1]->generate_tex(file);
            }
            else if (s == "/") {
                file << " \\frac{";
                c[0]->strip_paren()->generate_tex(file);
                file << "}{";
                c[1]->strip_paren()->generate_tex(file);
                file << "} ";
            }
            else if (s == "op") {
                c[0]->generate_tex(file);
                file << " \\otimes ";
                c[1]->generate_tex(file);
            }
            else if (s == "del") {
                file << " \\nabla_{";
                c[0]->generate_tex(file);
                file << "}{";
                c[1]->generate_tex(file);
                file << "} ";
            }
            else if (s=="^" || s=="_") {
                file << "{";
                c[0]->generate_tex(file);
                file << "}";
                file << s;
                file << "{";
                c[1]->generate_tex(file);
                file << "}";
            }
            else if (s == "lt") {
                c[0]->generate_tex(file);
                file << " < ";
                c[1]->generate_tex(file);
            }
            else if (s == "gt") {
                c[0]->generate_tex(file);
                file << " < ";
                c[1]->generate_tex(file);
            }
            else if (s == " le") {
                c[0]->generate_tex(file);
                file << " \\leq ";
                c[1]->generate_tex(file);
            }
            else if (s == "ge") {
                c[0]->generate_tex(file);
                file << " \\geq ";
                c[1]->generate_tex(file);
            }
            else if (s == "eq") {
                c[0]->generate_tex(file);
                file << " == ";
                c[1]->generate_tex(file);
            }
            else if (s == "ne") {
                c[0]->generate_tex(file);
                file << " \\neq ";
                c[1]->generate_tex(file);
            }
            else {
                c[0]->generate_tex(file);
                if (!(s=="+" || s=="-" || s=="*")) file << "\\ \\textsf{";
                file << s;
                if (!(s=="+" || s=="-" || s=="*")) file << "}\\ ";
                c[1]->generate_tex(file);
            }
        }
        else if (c.size() == 3) {
            // def      variable = paramlist -> exp
            // domain   lo, hi, D
            // ^_       x^y_z

            if (s == "def") {
                c[0]->generate_tex(file);
                file << " \\ = \\ ";
                c[1]->generate_tex(file);
                file << " \\rightarrow ";
                c[2]->generate_tex(file);
            }
            else if (s == "domain") {
                file << "\\Omega \\ = \\ [";
                c[0]->generate_tex(file);
                file << " , ";
                c[1]->generate_tex(file);
                file << "]^";
                c[2]->generate_tex(file);
            }
            else if (s == "^_") {
                c[0]->generate_tex(file);
                file << "^{";
                c[1]->generate_tex(file);
                file << "}_{";
                c[2]->generate_tex(file);
                file << "} ";
            }
            else if (s == "for") {
                file << " \\textsf{for}\\ ";
                c[0]->generate_tex(file);
                file << " \\in [";
                c[1]->generate_tex(file);
                file << ",";
                c[2]->generate_tex(file);
                file << "]";
            }
            else {
                std::cout << "Trying to generate unknown expression with 3 arguments?\n";
                this->print_tree(std::cout);
            }
        }
        else if (s == "sum") {
            file << " \\sum_{";
            c[0]->generate_tex(file);
            file << "=";
            c[1]->generate_tex(file);
            file << "}^{";
            c[2]->generate_tex(file);
            file << "} \\left(";
            c[3]->generate_tex(file);
            file << "\\right) ";
        }
        else {
            std::cout << "Trying to generate expression with more than 3 arguments?\n";
            this->print_tree(std::cout);
        }
    }

    const Exp* strip_paren() const {
        if (s == "()") return c[0];
        else return this;
    }

    void generate_cxx(std::ostream& file) const {
        // varlist, arglist, paramlist can be 1 or more entries
        if (s == "arglist" || s=="varlist" || s=="paramlist") {
            const_iterator iter = c.begin(); 
            while (iter != c.end()) {
                const Exp* e = *iter;
                
                if (s == "paramlist") {
                    Exp::Type t = e->type();
                    bool pass_by_const_ref = (! Exp::is_pod(t));
                    if (pass_by_const_ref) file << "const ";
                    file << Exp::cxx_type_name(t);
                    if (pass_by_const_ref) file << "&";
                    file << " ";
                }
                e->generate_cxx(file);
                ++iter;
                if (iter == c.end()) {
                    break;
                }
                else {
                    file << ", ";
                }
            }
        }
        else if (c.size() == 0) {
            if (s=="Let") {
                file << "#define WORLD_INSTANTIATE_STATIC_TEMPLATES\n";
                file << "#include <mra/mra.h>\n";
                file << "using namespace madness;\n";
            }
            else if (s == "In") {
                file << "\nvoid setup(World& world, int argc, char** argv) {\n";
                file << "    startup(world,argc,argv);\n";
                file << "    std::cout.precision(8);\n";
                file << "    FunctionDefaults<D>::set_cubic_cell(box_lo,box_hi);\n";
                file << "    FunctionDefaults<D>::set_k(k);\n";
                file << "    FunctionDefaults<D>::set_thresh(epsilon);\n";
                file << "    FunctionDefaults<D>::set_truncate_on_project(true);\n";
                file << "    FunctionDefaults<D>::set_project_randomize(true);\n";
                //file << "    FunctionDefaults<D>::set_truncate_mode(1);\n";
                file << "    real_grad = gradient_operator<double,D>(world);\n";
                file << "    complex_grad = gradient_operator<double_complex,D>(world);\n";
                file << "}\n";
                file << "\n";
                file << "int main(int argc, char**argv) {\n";
                file << "    initialize(argc,argv);\n";
                file << "    World world(MPI::COMM_WORLD);\n";
                file << "    setup(world, argc, argv);\n";
            }
            else if (s == "End") {
                file << "\n";
                file << "    real_grad.clear();\n";
                file << "    complex_grad.clear();\n";
                file << "    finalize();\n";
                file << "    return 0;\n";
                file << "}\n";
            }
            else if (s == "pi") {
                file << "constants::pi";
            }
            else if (s != "end") {
                file << s;
            }
        }
        else if (c.size() == 1) {
            // ()       (exp)
            // {}       {exp}
            // -        - exp
            // not      not exp
            // F        F variable
            // if, elif, while
            // greek    greek name
            // print    print arglist
            // plot     plot  exp
            // lap
            // trace

            if      (s == "()") file << "(";
            else if (s == "{}") file << "{";
            else if (s == "-") file << "-";
            else if (s == "F") {
                if (is_real(t)) file << "real_function(real_factory(world).f(";
                else            file << "complex_function(complex_factory(world).f(";
            }
            else if (s == "lap") file << "laplacian(";
            else if (s == "print") file << "print(";
            else if (s == "plot")  file << "plot(";
            else if (s == "if")  file << "if (";
            else if (s == "elif")  file << "else if (";
            else if (s == "while")  file << "while (";
            else if (s == "greek") ;
            else if (s == "not") file << "!";
            else if (s == "trace") file << "(";
            else file << s << " ";

            c[0]->generate_cxx(file);

            if (s == "{}") file << "}";
            if (s == "()" || s == "print" || s == "plot" || s == "if" || s=="elif" || s=="while" || s=="lap") file << ")";
            else if (s == "F") file << "))";
            else if (s == "()") file << ")";
            else if (s == "trace") file << ").trace()";
        }
        else if (c.size() == 2) {
            // decl     varlist in field
            // call     variable(arglist)
            // mathfun  funcname(exp)
            // braket   <a|b>
            // =, +, -, *, _, ^    exp op exp
            // comparison/logical ops
            // del
            // bsh

            if (s == "decl") {
            }
            else if (s == "call") {
                c[0]->generate_cxx(file);
                file << "(";
                c[1]->generate_cxx(file);
                file << ")";
            }
            else if (s == "mathfun") {
                c[0]->generate_cxx(file);
                file << "(";
                c[1]->generate_cxx(file);
                file << ")";
            }
            else if (s == "braket") {
                file << "inner(";
                c[0]->generate_cxx(file);
                file << ", ";
                c[1]->generate_cxx(file);
                file << ")";
            }
            else if (s == "=") {
                std::string cs = c[0]->str();
                if (cs == "greek") cs = c[0]->c[0]->str();
                // Declare variables at first assignment
                if (dectab.find(cs) == dectab.end()) {
                    dectab[cs] = true;
                    file << Exp::cxx_type_name(c[0]->type()) << " ";
                }
                c[0]->generate_cxx(file);
                file << " = ";
                c[1]->generate_cxx(file);
            }
            else if (s == "/" && is_integer(c[0]->type()) && is_integer(c[1]->type())) {
                // Probably want 1/2=0.5
                file << "double(";
                c[0]->generate_cxx(file);
                file << ") / ";
                c[1]->generate_cxx(file);
            }
            else if (s == "^") {
                if (c[1]->str() == "2") {
                    file << "(";
                    c[0]->generate_cxx(file);
                    file << "*";
                    c[0]->generate_cxx(file);
                    file << ")";
                }
                else if (c[1]->str()=="-" && c[1]->c.size()==1 && c[1]->c[0]->str()=="1") {
                    file << "(1.0/";
                    c[0]->generate_cxx(file);
                    file << ")";
                }
                else {
                    file << "std::pow(";
                    c[0]->generate_cxx(file);
                    file << ", ";
                    c[1]->generate_cxx(file);
                    file << ")";
                }
            }
            else if (s == "_") {
                c[0]->generate_cxx(file);
                file << "[";
                c[1]->generate_cxx(file);                
                file << "]";
            }            
            else if (s == "and") {
                c[0]->generate_cxx(file);
                file << " && ";
                c[1]->generate_cxx(file);
            }
            else if (s == "or") {
                c[0]->generate_cxx(file);
                file << " || ";
                c[1]->generate_cxx(file);
            }
            else if (s == "lt") {
                c[0]->generate_cxx(file);
                file << " < ";
                c[1]->generate_cxx(file);
            }
            else if (s == "gt") {
                c[0]->generate_cxx(file);
                file << " < ";
                c[1]->generate_cxx(file);
            }
            else if (s == "le") {
                c[0]->generate_cxx(file);
                file << " <= ";
                c[1]->generate_cxx(file);
            }
            else if (s == "ge") {
                c[0]->generate_cxx(file);
                file << " >= ";
                c[1]->generate_cxx(file);
            }
            else if (s == "eq") {
                c[0]->generate_cxx(file);
                file << " == ";
                c[1]->generate_cxx(file);
            }
            else if (s == "ne") {
                c[0]->generate_cxx(file);
                file << " != ";
                c[1]->generate_cxx(file);
            }
            else if (s == "del") {
                if (is_real(t)) {
                    file << "(*real_grad[";
                    c[0]->generate_cxx(file);
                    file << "])(";
                    c[1]->generate_cxx(file);
                    file << ")";
                }
            }
            else if (s == "bsh") {
                file << "BSH(sqrt(";
                c[0]->generate_cxx(file);
                file << "))((";
                c[1]->generate_cxx(file);
                file << ").truncate()).truncate()";
                if (c[0]->str() == "0.0") file <<".scale(-1.0)"; // Since mapped del^-2 not -del^-2
            }
            else {
                c[0]->generate_cxx(file);
                file << " " << s << " ";
                c[1]->generate_cxx(file);
            }
        }
        else if (c.size() == 3) {
            // def      variable = paramlist -> exp
            // domain   lo, hi, D
            // ^_       x^y_z

            if (s == "def") {
                file << "\n" << Exp::cxx_type_name(c[0]->type()) << " ";
                c[0]->generate_cxx(file);
                file << "(";
                c[1]->generate_cxx(file);
                file << ") {\n    return ";
                c[2]->generate_cxx(file);
                file << ";\n}";
            }
            else if (s == "domain") {
                file << "\nconst double box_lo = ";
                c[0]->generate_cxx(file);
                file << ", box_hi = ";
                c[1]->generate_cxx(file);
                file << ";\n";
                file << "const int D = ";
                c[2]->generate_cxx(file);
                file << ";\n";
                file << "typedef Function<double,D> real_function;\n";
                file << "typedef Function<double_complex,D> complex_function;\n";
                file << "typedef FunctionFactory<double,D> real_factory;\n";
                file << "typedef FunctionFactory<double_complex,D> complex_factory;\n";
                file << "typedef Vector<double,D> coord;\n\n";
                file << "std::vector< std::shared_ptr < Derivative<double,D> > > real_grad;\n";
                file << "std::vector< std::shared_ptr < Derivative<double_complex,D> > > complex_grad;\n";
                file << "template <typename functionT>\n";
                file << "void plot(const functionT& f, const char* filename = \"plot.dx\") {\n";
                file << "    plotdx(f, filename, FunctionDefaults<D>::get_cell(), std::vector<long>(D,101));\n";
                file << "}\n";
                file << "template <typename functionT> functionT laplacian(const functionT& f) {throw \"not yet\";}\n";
                file << "template <typename T> double norm(const T& t) {return t.normf();}\n";
                file << "template <typename T> double normf(const T& t) {return t.normf();}\n";
                file << "template <typename T> double norm(const Function<T,D>& t) {return t.norm2();}\n";
                file << "#define BSH(mu) BSHOperator<D>(world, mu, 1e-4, epsilon)";

            }
            else if (s == "^_") {
                if (c[1]->str() == "2") {
                    file << "(";
                    c[0]->generate_cxx(file);
                    file << "[";
                    c[2]->generate_cxx(file);
                    file << "]*";
                    c[0]->generate_cxx(file);
                    file << "[";
                    c[2]->generate_cxx(file);
                    file << "]";
                    file << ")";
                }
                else {
                    file << "std::pow(";
                    c[0]->generate_cxx(file);
                    file << "[";
                    c[2]->generate_cxx(file);
                    file << "], ";
                    c[1]->generate_cxx(file);
                    file << ") ";
                }
            }
            else if (s == "for") {
                file << "for (int ";
                c[0]->generate_cxx(file);
                file << "=";
                c[1]->generate_cxx(file);
                file << "; ";
                c[0]->generate_cxx(file);
                file << "<=";
                c[2]->generate_cxx(file);
                file << "; ";
                c[0]->generate_cxx(file);
                file << "++) {";
            }
            else {
                std::cout << "Trying to generate unknown expression with 3 arguments?\n";
                this->print_tree(std::cout);
            }
        }
        else if (s == "sum") {
            file << "for (int ";
            c[0]->generate_cxx(file);
            file << "=";
            c[1]->generate_cxx(file);
            file << "; ";
            c[0]->generate_cxx(file);
            file << "<=";
            c[2]->generate_cxx(file);
            file << "; ";
            c[0]->generate_cxx(file);
            file << "++) {";
            c[3]->generate_cxx(file);
            file << " ;}";
        }
        else {
            std::cout << "Trying to generate expression with more than 3 arguments?\n";
            this->print_tree(std::cout);
        }
    }

    const std::string& str() const {
        return s;
    }

    Type type() const {
        return t;
    }

    iterator begin() {
        return c.begin();
    }

    iterator end() {
        return c.end();
    }

    const_iterator begin() const {
        return c.begin();
    }

    const_iterator end() const {
        return c.end();
    }

    // Type of MRA function inferred from C++ function return type
    static Type funtype(Type t) {
        if      (t == REAL) return REALFUN;
        else if (t == COMPLEX) return COMPLEXFUN;
        else {
            std::cout << "tying to make a function from " << type_name(t) << std::endl;
            return UNKNOWN;
        }
    }

    // Type of result from subscripting/indexing operation t_i == t[i]
    static Type subscrtype(Type t) {
        if (t==COORD || t==REALTEN) return REAL;
        else if (t==COMPLEXTEN) return COMPLEX;
        else return UNKNOWN;
    }

    // Type of result from braket <f|g>
    static Type brakettype(Type f, Type g) {
        if (f==REALFUN) {
            if (g==REALFUN) return REAL;
            else if (g==COMPLEXFUN) return COMPLEX;
        }
        else if (f==COMPLEXFUN && (g==REALFUN || g==COMPLEXFUN)) 
            return COMPLEX;
        
        std::cout << "computing braket with non-function types? <" << type_name(f) << "|" << type_name(g) << ">\n";
        return UNKNOWN;
    }

    static bool is_real(Type t) {
        return (t==INTEGER || t==REAL || t==REALFUN || t==REALTEN);
    }

    static bool is_complex(Type t) {
        return (t==COMPLEX || t==COMPLEXFUN || t==COMPLEXTEN);
    }

    static bool is_string(Type t) {
        return (t==STRING);
    }

    static bool is_tensor(Type t) {
        return (t==REALTEN || t==COMPLEXTEN);
    }

    static bool is_integer(Type t) {
        return (t==INTEGER);
    }

    static bool is_number(Type t) {
        return (t==INTEGER || t==REAL);
    }

    static bool is_unknown(Type t) {
        return (t==UNKNOWN);
    }

    static bool is_pod(Type t) {
        return (t==INTEGER || t==REAL || t==COMPLEX);
    }

    static bool is_mrafun(Type t) {
        return (t==REALFUN || t==COMPLEXFUN);
    }

    // Try to guess type of result from unaryop(op,t)
    static Type unaryoptype(const std::string& op, Type t) {
        if (op=="norm" || op=="normf") return REAL;
        return t;
    }

    // Try to guess type of result from binaryop(op,l,r)
    static Type binaryoptype(const std::string& op, Type l, Type r) {
        if (op == "mathfun") return unaryoptype(op, r);
        if (op == "norm" || op == "normf") return REAL;
        if (is_unknown(l) || is_unknown(r)) return UNKNOWN;
        if (is_string(l) || is_string(r)) return STRING;

        // Can it be an integer?
        if (is_integer(l) && is_integer(r)) {
            if (op=="+" || op=="-" || op=="*") return INTEGER;
        }
        
        // Choice is now (complex|real) (number|tensor|function)
        if (is_complex(l) || is_complex(r)) {
            if (is_tensor(l) || is_tensor(r)) return COMPLEXTEN;
            if (is_mrafun(l) || is_mrafun(r)) return COMPLEXFUN;
            return COMPLEX;
        }
        else {
            if (is_tensor(l) || is_tensor(r)) return REALTEN;
            if (is_mrafun(l) || is_mrafun(r)) return REALFUN;
            return REAL;
        }
    }

    void set_type(Type t) {
        this->t = t;
    }

    ~Exp() {
        for (iterator iter = c.begin(); 
             iter != c.end(); 
             ++iter) delete *iter;
    }

    Exp* operator[](int i) const {
        return c[i];
    }
    
};

inline std::ostream& operator<<(std::ostream& s, const Exp& e) {
    s << Exp::type_name(e.type()) << " " << e.str(); 
    return s;
}

// Conducting the whole scanning and parsing of mra.
class mra_driver {
public:
    int scopedepth;             // Nesting depth of if/while to disambiguate end keyword
    int tmpvarcnt;              // Counts temporary variables
    bool use_k_default;         // False if user specified k
    bool use_eps_default;       // False if user specified eps

    std::vector<Exp*>::iterator cur; // When traversing entire tree this is the current statement
    std::map<std::string, Exp::Type> symtab;   // maps symbols to types
    std::vector<Exp*> program;                 // The full AST


    mra_driver ();
    virtual ~mra_driver ();


    // Handling the scanner.
    void scan_begin ();
    void scan_end ();
    bool trace_scanning;

    std::string tmpvar() {
        std::ostringstream o;
        o << "tmp" << tmpvarcnt++;
        return o.str();
    }

    // Run the parser.  Return 0 on success.
    int parse (const std::string& f);
    std::string file;
    bool trace_parsing;

    // Error handling.
    void error (const yy::location& l, const std::string& m);
    void error (const std::string& m);

    // Adds a new statement
    Exp* add(Exp* statement) {
        program.push_back(statement);
        return statement;
    }

    // Inserts a new statement *before* current one and maintains current iterator
    void insert_before(Exp* e) {
        cur = program.insert(cur, e);
        cur++;
    }

    // Inserts a new statement *after* current one and maintains current iterator
    void insert_after(Exp* e) {
        cur = program.insert(cur+1, e);
        cur--;
    }

    // Rewrites statement to replace sum with for loops prior to generating C++
    void rewrite_sum(Exp* ee) {
        Exp& e = *ee;
        const std::string& s = e.str();
        const Exp::Type t = e.type();
        if (s == "sum") {
            Exp* ind = e[0];
            Exp* lo  = e[1];
            Exp* hi  = e[2];
            Exp* exp = e[3];

            //if (Exp::is_pod(t)) {
            Exp* var = new Exp(tmpvar(),t);
            insert_sym("temporary", var);
            insert_before(new Exp("=", var, new Exp("0",Exp::INTEGER)));
            insert_before(new Exp("for",ind,lo,hi));

            rewrite_sum(exp);

            insert_before(new Exp("=",new Exp(*var), 
                                  new Exp("+",new Exp(*var),exp,t),t));
            insert_before(new Exp("end"));
            e = *var;
            //}
            
        }
        else {
            for (Exp::iterator it=e.begin(); it!=e.end(); ++it) {
                rewrite_sum(*it);
            }
        }
    }

    Exp::Type lookup_type(const std::string& s) {
        std::map<std::string,Exp::Type>::const_iterator it = symtab.find(s);
        if (it != symtab.end()) {
            return it->second;
        }
        else if (s[0] == 'x') { /* Implicit typing of coordinates x */
            std::cout << "implict typing of " << s << " as Rd (coordinate)\n";
            symtab[s] = Exp::COORD;
            return Exp::COORD;
        }
        else {
            return Exp::UNKNOWN;
        }
    }

    Exp::Type lookup_type(const Exp* e) {
        return lookup_type(e->str());
    }

    void insert_sym(const std::string& how, const std::string& s, Exp::Type t) {
        Exp::Type oldt = lookup_type(s);
        if (oldt != t && oldt != Exp::UNKNOWN) {
            std::cout << how << " type of " << s << " as " << Exp::type_name(t) << " overriding previous definition of " << Exp::type_name(oldt) << std::endl;
        }
        else {
            std::cout << how << " type of " << s << " as " << Exp::type_name(t) << std::endl;
        }
        symtab[s] = t;
    }

    void insert_sym(const std::string& how, const Exp* e) {
        if (e->str() == "greek") {
            insert_sym(how, (*(e->begin()))->str(), e->type()); // gack!
        }
        else {
            insert_sym(how, e->str(), e->type());
        }
    }

    void addsyms(const Exp* e, const std::string& type) {
        Exp::Type t;

        if      (type == "Z")  t = Exp::INTEGER;
        else if (type == "R")  t = Exp::REAL;
        else if (type == "C")  t = Exp::COMPLEX;
        else if (type == "TR") t = Exp::REALTEN;
        else if (type == "TC") t = Exp::COMPLEXTEN;
        else if (type == "FR") t = Exp::REALFUN;
        else if (type == "FC") t = Exp::COMPLEXFUN;
        else if (type == "S")  t = Exp::STRING;
        else {
            std::cout << "unknown type " << type << " in declaration" << std::endl;
            t = Exp::UNKNOWN;
        }
        
        if (e->str() == "varlist") {
            for (Exp::const_iterator iter=e->begin(); iter!=e->end(); ++iter) {
                insert_sym("declaring",(*iter)->str(), t);
            }
        }
        else {
            std::cout << "!! was expecting varlist in declaration\n";
        }
    }

    void print_tree(std::ostream& file) const {
        for (std::vector<Exp*>::const_iterator iter = program.begin(); iter != program.end();  ++iter) 
            (*iter)->print_tree(file);
    }

    void regenerate(std::ostream& file) const {
        int indent = 1;
        for (std::vector<Exp*>::const_iterator iter = program.begin(); iter != program.end();  ++iter) {
            const Exp* e = *iter;
            const std::string& s = e->str();

            int ind = indent;
            if (s=="Let" || s=="In" || s=="End") {
                ind = 0;
            }
            else if (s=="if" || s=="while" || s=="for") {
                indent++;
            }
            else if (s=="elif" || s=="else") {
                ind--;
            }
            else if (s == "end") {
                ind--;
                indent--;
            }

            for (int i=0; i<4*ind; i++) file << " ";
            
            (*iter)->regenerate(file);

            file << std::endl;
        }
    }

    void generate_tex(std::ostream& file) const {
        int indent = 1;
        for (std::vector<Exp*>::const_iterator iter = program.begin(); iter != program.end();  ++iter) {
            const Exp* e = *iter;
            const std::string& s = e->str();

            int ind = indent;
            if (s=="Let" || s=="In" || s=="End") {
                ind = 0;
            }
            else if (s=="if" || s=="while" || s=="for") {
                indent++;
            }
            else if (s=="elif" || s=="else") {
                ind--;
            }
            else if (s == "end") {
                ind--;
                indent--;
            }

            file << "& ";
            for (int i=0; i<4*ind; i++) file << "\\ ";
            
            (*iter)->generate_tex(file);

            file << "\\\\" << std::endl;
        }
    }

    void generate_cxx(std::ostream& file) {

        // note class global variable cur
        for (cur = program.begin(); cur != program.end();  ++cur) 
            rewrite_sum(*cur);

        int indent = 0;
        for (std::vector<Exp*>::const_iterator iter = program.begin(); iter != program.end();  ++iter) {
            const Exp* e = *iter;
            const std::string& s = e->str();

            int ind = indent;
            if (s=="Let" || s=="In" || s=="End") {
                ind = 0;
                if (s=="In") indent++;
            }
            else if (s=="if" || s=="while" || s=="for") {
                indent++;
            }
            else if (s=="elif" || s=="else") {
                ind--;
            }
            else if (s == "end") {
                ind--;
                indent--;
            }

            if (s == "=") {
                const std::string& cs = (*(e->begin()))->str();
                if (cs == "k") use_k_default = false;
                if (cs == "greek" && (*(*(e->begin()))->begin())->str()=="epsilon") 
                    use_eps_default = false;
            }

            if (s == "In") {
                if (use_k_default) {
                    file << "const int k = 6;\n";
                    insert_sym("defaulting", new Exp("k",Exp::INTEGER));
                }
                if (use_eps_default) {
                    file << "const double epsilon = 1e-4;\n";
                    insert_sym("defaulting", new Exp("epsilon",Exp::REAL));
                }
            }                


            for (int i=0; i<4*ind; i++) file << " ";
            
            if (s=="elif" || s=="else") {
                file << "}\n";
                for (int i=0; i<4*ind; i++) file << " ";
            }
            
            (*iter)->generate_cxx(file);

            if (s=="if" || s=="elif" || s=="else" || s=="while") 
                file << " {";
            else if (s=="end")
                file << "}";
            else if (!(s=="Let" || s=="In" || s=="End" || s=="def" || s=="domain" || s=="sum" || s=="for"))
                file << ";";

            file << std::endl;
        }
    }

};

#endif // ! MRA_DRIVER_HH
