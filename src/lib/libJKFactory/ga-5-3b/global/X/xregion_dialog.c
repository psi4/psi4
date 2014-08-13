#if HAVE_CONFIG_H
#   include "config.h"
#endif

/* This module contains functions for creating and destroying a simple
   dialog box.  
*/

#include "xregion.h"

void dialog_box(char *dlg_str)
{
  Display *disp;

  disp = XtDisplay(overview_shell);

  /* Create top level shell widget */
  dlg_top = XtVaAppCreateShell("xregion","XRegion",
                     applicationShellWidgetClass,disp,
                     NULL);

  /* Create form widget to hold everything else */

  dlg_form = XtVaCreateManagedWidget("dialogform", formWidgetClass,
                                     dlg_top,
                                     XtNbackground, DEFAULT_BG,
                                     XtNforeground, DEFAULT_FG,
                                     NULL);

  dlg_label = XtVaCreateManagedWidget("dialoglabel", labelWidgetClass,
                                     dlg_form,
                                     XtNbackground, DEFAULT_BG,
                                     XtNforeground, DEFAULT_FG,
                                     XtNlabel, dlg_str,
                                     XtNvertDistance, 5,
                                     XtNhorizDistance, 5,
                                     XtNborderWidth, 0,
                                     NULL);

  dlg_btn = XtVaCreateManagedWidget("dialogbutton", commandWidgetClass,
                                     dlg_form,
                                     XtNbackground, DEFAULT_BG,
                                     XtNforeground, DEFAULT_FG,
                                     XtNlabel, "Dismiss",
                                     XtNfromVert, dlg_label,
                                     XtNvertDistance, 5,
                                     XtNhorizDistance, 30,
                                     NULL);
  XtAddCallback(dlg_btn, XtNcallback, dismiss_dialog, NULL);
 
  XtRealizeWidget(dlg_top);
}

/* JJU: void dismiss_dialog(Widget w, caddr_t data, XEvent *event) */
void dismiss_dialog(Widget w, XtPointer data, XtPointer event)
{
  /* remove dialog box */
  XtDestroyWidget(dlg_btn);
  XtDestroyWidget(dlg_label);
  XtDestroyWidget(dlg_form);
  XtDestroyWidget(dlg_top);
}

