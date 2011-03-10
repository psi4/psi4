namespace opt {

/**
  \ingroup optking
  \brief   Exception class for problems with internal coordinates.
  */

/* If they relate to values of the coordinates and derivatives, then try new
   coordinates ; if it looks like user error in definition, then
   quit right away.  */
class INTCO_EXCEPT {
  private:
   const char * message;
   bool try_other_intcos;

  public:
    INTCO_EXCEPT(const char * m) {
      message = m;
      try_other_intcos = false;
    }

    INTCO_EXCEPT(const char * m, bool t) {
      message = m;
      try_other_intcos = t;
    }

    ~INTCO_EXCEPT() {};

    bool try_again() { return try_other_intcos; }
    const char *g_message(void) { return message; }
};

class BAD_STEP_EXCEPT {
  private:
   const char * message;

  public:
    BAD_STEP_EXCEPT(const char * m) {
      message = m;
    }

    BAD_STEP_EXCEPT(const char * m, bool t) {
      message = m;
    }

    ~BAD_STEP_EXCEPT() {};

    const char *g_message(void) { return message; }
};


}

