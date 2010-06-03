/* $Id: plotQt.c 7523 2005-12-09 19:34:24Z kb $

Copyright (C) 2000  The PARI group.

This file is part of the PARI/GP package.

PARI/GP is free software; you can redistribute it and/or modify it under the
terms of the GNU General Public License as published by the Free Software
Foundation. It is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY WHATSOEVER.

Check the License for details. You should have received a copy of it, along
with the package; see the file 'COPYING'. If not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA. */
/////////////////////////////////////////////////////////////////////////////
//
//  High resolution plot using Trolltech's Qt library
//
//  You may possibly want to use this file with a "Qt Free Edition"
//  which is distributed under the terms of the Q PUBLIC LICENSE (QPL),
//  or with a "Qt/Embedded Free Edition" which is
//  distributed under the terms of the GNU General Public License (GPL).
//  Please check http://www.trolltech.com for details.
//
//                            ---Nils-Peter Skoruppa (www.countnumber.de)
/////////////////////////////////////////////////////////////////////////////
extern "C" {
#include "pari.h"
#undef grem
#include "rect.h"
}

#ifdef __QPE__
#include <qpeapplication.h>
#else
#include <qapplication.h>
#endif
#include <qwidget.h>
#include <qpainter.h>
#include <qarray.h>
#include <qpoint.h>
#include <qrect.h>
#include <qcolor.h>
#include <qpixmap.h>
#include <qimage.h>


class Plotter: public QWidget {

#ifdef __FANCY_WIN__
     Q_OBJECT

signals:
    void clicked();

protected:
    void mouseReleaseEvent( QMouseEvent*);
#endif

public:
    Plotter( long *w, long *x, long *y, long lw,
	     QWidget* parent = 0, const char* name = 0, WFlags fl = 0);
    void save( const QString& s = *plotFile + ".xpm",//QString("pariplot.xpm"),
	       const QString& f = QString( "XPM"));

protected:
    void paintEvent( QPaintEvent *);
    void resizeEvent ( QResizeEvent *);
#ifndef __FANCY_WIN__
    void keyPressEvent( QKeyEvent *);
#endif

private:
    long *w;                           // map into rectgraph indexes
    long *x;                           // x, y: array of x,y-coorinates of the
    long *y;                           //   top left corners of the rectwindows
    long lw;                           // lw: number of rectwindows
    QColor color[MAX_COLORS];
    QFont font;
    static QString *plotFile;
    void draw(QPainter *p);

// public:
//     static void setPlotFile( const char *);
};


QString *Plotter::plotFile = new QString( "pariplot");


Plotter::Plotter( long *w, long *x, long *y, long lw,
		  QWidget* parent,  const char* name, WFlags fl)
    : QWidget( parent, name, fl), font( "lucida", 9) {

    this->w=w; this->x=x; this->y=y; this->lw=lw;
#ifndef __FANCY_WIN__
    this->resize( pari_plot.width, pari_plot.height);
    this->setCaption( "Pari QtPlot");
#endif
    this->setBackgroundColor( white);
    this->setFont( font);
    // defaults as in plotX.c (cf. rgb.txt in a X11 system)
    color[0]         = white;
    color[BLACK]     = black;
    color[BLUE]      = blue;
    color[VIOLET]    = QColor( 208, 32, 144) ;
    color[RED]       = red;
    color[GREEN]     = green;
    color[GREY]      = gray;
    color[GAINSBORO] = QColor( 220, 220, 220);
}

// void Plotter::setPlotFile( const char *s) {

//     delete Plotter::plotFile;
//     Plotter::plotFile = new QString( s);
// }

struct data_qt
{
  QPainter *p;
  QColor *color;
};

static void SetForeground(void *data, long col)
{
   struct data_qt *d = (struct data_qt *) data;
   d->p->setPen(d->color[col]);
}

static void DrawPoint(void *data, long x, long y)
{
   struct data_qt *d = (struct data_qt *) data;
   d->p->drawPoint(x, y);
}

static void DrawLine(void *data, long x1, long y1, long x2, long y2)
{
   struct data_qt *d = (struct data_qt *) data;
   d->p->drawLine(x1, y1, x2, y2);
}

static void DrawRectangle(void *data, long x, long y, long w, long h)
{
   struct data_qt *d = (struct data_qt *) data;
   d->p->drawRect(x, y, w, h);
}

static void DrawPoints(void *data, long nb, struct plot_points *p)
{
   struct data_qt *d = (struct data_qt *) data;
   QPointArray xp=QPointArray(nb);
   long i;
   for (i=0;i<nb;i++)
     xp.setPoint(i,p[i].x, p[i].y);
   d->p->drawPoints(xp);
}

static void DrawLines(void *data, long nb, struct plot_points *p)
{
   struct data_qt *d = (struct data_qt *) data;
   QPointArray xp=QPointArray(nb);
   long i;
   for (i=0;i<nb;i++)
     xp.setPoint(i,p[i].x, p[i].y);
   d->p->drawPolyline(xp);
}

static void DrawString(void *data, long x, long y, char *text, long numtext)
{
  struct data_qt *d = (struct data_qt *) data;
  d->p->drawText(x, y, QString(text), numtext);
}

void Plotter::draw(QPainter *p){
  struct plot_eng plotQt;
  struct data_qt d;
  d.p= p;
  d.color=color;
  plotQt.sc=&SetForeground;
  plotQt.pt=&DrawPoint;
  plotQt.ln=&DrawLine;
  plotQt.bx=&DrawRectangle;
  plotQt.mp=&DrawPoints;
  plotQt.ml=&DrawLines;
  plotQt.st=&DrawString;
  plotQt.pl=&pari_plot;
  double xs = double(this->width()) / pari_plot.width,
         ys = double(this->height()) / pari_plot.height;
  gen_rectdraw0(&plotQt, (void *)&d, this->w, this->x, this->y,this->lw,xs,ys);
}

void Plotter::save( const QString& s, const QString& f){

    QPixmap pm( this->width(), this->height());
    QPainter p;

    p.begin( &pm, this);
    p.fillRect( 0, 0, pm.width(), pm.height(), white);
    this->draw(&p);
    p.end();

    // supported formats in qt2: BMP, JPEG, PNG, PNM, XBM, XPM ; PNG is broken
    pm.save( s, f);
}

void Plotter::paintEvent( QPaintEvent *) {

    QPainter p;
    p.begin( this);
    this->draw(&p);
    p.end();
}

void Plotter::resizeEvent( QResizeEvent *) { }

#ifndef __FANCY_WIN__
void Plotter::keyPressEvent( QKeyEvent *e) {

    switch ( tolower( e->ascii())) {
        case 's':
	    save();
	    this->setCaption( "Pari QtPlot: " + *plotFile);
            break;
    }
}
#endif


#ifdef __FANCY_WIN__
void Plotter::mouseReleaseEvent( QMouseEvent*) {

    emit clicked();
}
#endif



#ifdef __FANCY_WIN__
//
// The envelope for an instance of plotter
//


#include <qmainwindow.h>
#include <qpopupmenu.h>
#include <qmenubar.h>
#include <qtoolbar.h>
#include <qaction.h>
#include <qfiledialog.h>
#include <qmessagebox.h> 
#include <qfile.h> 
#include <qstatusbar.h>
#include <qimage.h>
#include <qstrlist.h>
#include <qlabel.h>
#include <qspinbox.h>
#include <qlayout.h>


/* XPM */
static const char * const fullscreen_xpm[] = {
"14 14 2 1",
" 	c None",
".	c #000000",
"..............",
".     ..     .",
".     ..     .",
".    ....    .",
".     ..     .",
".  .  ..  .  .",
"..............",
"..............",
".  .  ..  .  .",
".     ..     .",
".    ....    .",
".     ..     .",
".     ..     .",
".............."}; 


class SaveAsDialog: public
#ifdef __QPE__
//QDialog
#else
QFileDialog
#endif
{

    Q_OBJECT

public:
    SaveAsDialog( const QString & c = QString::null,
		  const QString & s = QString::null, int w = 0, int h = 0,
		  QWidget *parent = 0, const char *name = 0, WFlags f = 0);
    ~SaveAsDialog();
#ifdef __QPE__
    QString selectedFile() { return nameW->text();}
#endif
    int picWidth() { return widthW->value();}
    int picHeight() { return heightW->value();}

private:
    QLineEdit *nameW;
    QSpinBox *widthW, *heightW;

};


SaveAsDialog::SaveAsDialog( const QString & c, const QString & s, int w, int h,
			    QWidget *parent, const char *name, WFlags f)
#ifdef __QPE__
    // simplistic dialog in case of QPE ( fancy alternative: class FileSelector)

    : QDialog( parent, name, TRUE, f) {

    if( c) this->setCaption( c);
    nameW = new QLineEdit( this);
    if( s) nameW->setText( s);
    widthW = new QSpinBox( 1, 65536, 1, this);
    if( w > 0) widthW->setValue( w);
    heightW = new QSpinBox( 1, 65536, 1, this);
    if( h > 0) heightW->setValue( h);

    QVBoxLayout *top = new QVBoxLayout( this, 10);
    QGridLayout *contents = new QGridLayout( 3, 2);

    top->addLayout( contents);

    QLabel *l;
    l = new QLabel( nameW, "Name : ", this);
    l->setAlignment( AlignRight | AlignVCenter);
    contents->addWidget( l, 0, 0);
    contents->addWidget( nameW, 0, 1);
    l = new QLabel( widthW, "Width : ", this);
    l->setAlignment( AlignRight | AlignVCenter);
    contents->addWidget( l, 1, 0);
    contents->addWidget( widthW, 1, 1);
    l = new QLabel( heightW, "Height : ", this);
    l->setAlignment( AlignRight | AlignVCenter);
    contents->addWidget( l, 2, 0);
    contents->addWidget( heightW, 2, 1);

    top->activate();
    this->resize( 160, this->height()); // hack!!!
#else
    : QFileDialog( parent, name, TRUE) {

    if( c) this->setFilters( c);
    if( s) this->setSelection( s);

    QLabel *l;
    QWidget *wt = new QWidget( this);
    QHBoxLayout *spinBoxes = new QHBoxLayout( wt, 5);
    widthW = new QSpinBox( 1, 65536, 1, wt);
    l  = new QLabel( widthW, "&width ", wt);
    spinBoxes->addWidget( l);
    spinBoxes->addWidget( widthW);
    if( w > 0) widthW->setValue( w);
    heightW = new QSpinBox( 1, 65536, 1, wt);
    spinBoxes->addSpacing(10);
    l  = new QLabel( heightW, "&height ", wt);
    l->setAlignment( AlignRight | AlignVCenter);
    spinBoxes->addWidget( l);
    spinBoxes->addWidget( heightW);
    if( h > 0) heightW->setValue( h);
    l = new QLabel( "Resolution:", this);
    QFileDialog::addWidgets( l, wt, 0);
#endif
}


SaveAsDialog::~SaveAsDialog() {
}



class PlotWindow: public QMainWindow {

     Q_OBJECT

public:
    PlotWindow( long *w, long *x, long *y, long lw,
		QWidget* parent = 0, const char* name = 0, WFlags fl = 0);
    ~PlotWindow();

#ifndef __QPE__
protected:
    void resizeEvent( QResizeEvent *);
#endif

private slots:
    void fullScreen();
    void normalView();
    void save();
    void save( int);

private:
    static const QStrList file_formats;
    Plotter *plr;
    QString saveFileName;
    int saveFileFormat;
#ifndef __QPE__
    QLabel *res;
#endif
};


const QStrList PlotWindow::file_formats = QImage::outputFormats();


PlotWindow::PlotWindow( long *w, long *x, long *y, long lw,
			QWidget* parent, const char* name, WFlags fl)
    : QMainWindow( parent, name, fl),
      saveFileName( "pariplot"), saveFileFormat( 0) {

    this->setCaption( "Pari QtPlot");

#ifdef __QPE__
    QToolBar *toolBar = new QToolBar( this);
    QMenuBar *menuBar = new QMenuBar( toolBar);
    toolBar->setHorizontalStretchable( TRUE);
    this->setToolBarsMovable( FALSE);
#else
    QMenuBar *menuBar = this->menuBar();
    menuBar->setFrameStyle( QFrame::Panel | QFrame::Raised);
#endif

    // Setting up the File and View menu
    QPopupMenu *format = new QPopupMenu( this);
    for( uint i = 0; i < file_formats.count(); i++) {
	format->insertItem( QString( QStrList(file_formats).at(i)) + " ...",
			    this, SLOT( save( int)), 0, i);
	if( 0 == QString( QStrList(file_formats).at(i)).compare( "PNG"))
	    format->setItemEnabled( i, FALSE); // PNG seems to be broken
    }
    QPopupMenu *file = new QPopupMenu( this);
    CHECK_PTR( file );
    file->insertItem( "&Save", this, SLOT( save()), CTRL+Key_S);
    file->insertItem( "Save &as", format);
    file->insertSeparator();
    file->insertItem( "&Close", this, SLOT(close()), CTRL+Key_C);
    menuBar->insertItem( "&File", file );

#ifndef __QPE__
    QPopupMenu *view = new QPopupMenu( this);
    menuBar->insertItem( "&View", view);
#endif
    // Setting up the Fullscreen action
    QAction *a = new QAction( "use full screen",
			      QPixmap( (const char ** )fullscreen_xpm),
			      "&Fullscreen", CTRL+Key_F, this);
    connect( a, SIGNAL( activated()), this, SLOT( fullScreen()));
#ifdef __QPE__
    a->addTo( toolBar);
#else
    a->addTo( view);
#endif

    // Setting up an instance of plotter
    plr = new Plotter( w, x, y, lw, this);
    connect( plr, SIGNAL(clicked()), this, SLOT( normalView()));
    this->setCentralWidget( plr);

#ifndef __QPE__
    this->resize( pari_plot.width,
		  pari_plot.height + 25);
    res = new QLabel( statusBar());
    statusBar()->addWidget( res);
#endif
}


PlotWindow::~PlotWindow() {
}


#ifndef __QPE__
void PlotWindow::resizeEvent( QResizeEvent *e) {
    
    QMainWindow::resizeEvent( e);
    res->setText( QString( "Resolution: ") +
		  QString::number( plr->width()) + "x" +
		  QString::number( plr->height()));
    res->setFixedSize( res->sizeHint());
}
#endif


void PlotWindow::fullScreen() {

    if ( plr->parentWidget()) {
	plr->reparent( 0, WStyle_Tool | WStyle_Customize | WStyle_StaysOnTop,
		      QPoint( 0, 0), FALSE);
	plr->resize( qApp->desktop()->width(), qApp->desktop()->height());
	plr->show();
    }
}


void PlotWindow::normalView() {

    if ( !plr->parentWidget()) {
	plr->reparent( this, 0, QPoint(0,0), FALSE);
	this->setCentralWidget( plr);
	plr->show();
    }
}


void PlotWindow::save() {

    QString ff = QString( QStrList(file_formats).at( saveFileFormat));
    QString fn = saveFileName + "." + ff.lower();
    plr->save( fn, ff);
    this->setCaption( QString( "Pari QtPlot:") + fn);
#ifndef __QPE__
    statusBar()->message( QString( "File %1 saved" ).arg( fn), 2000 );
#endif
}


void PlotWindow::save( int id) {

    QString ff( QStrList(file_formats).at( id));
#ifdef __QPE__
    QString s( "Save as");
#else
    QString s( ff + " (*." + ff.lower() +");;All (*)");
#endif
    SaveAsDialog d( s, saveFileName + "." + ff.lower(),
		    plr->width(), plr->height());
    if( QDialog::Rejected == d.exec()) return;
    QString fn = d.selectedFile();
    if ( !fn.isEmpty()) {
	if( QFile( fn).exists() &&
	    QMessageBox::warning(
		this, this->caption(),
		QString( "A file named\n\"") + fn
		+ QString( "\"\nalready exists\n"
			   "Should this file be overwritten ?\n\n"),
		"&Overwrite", "&Cancel"))  return;
	saveFileName = fn;
	int p;
	if( (p = saveFileName.findRev( "." + ff, -1, FALSE)) >=0)
	    saveFileName.truncate( p);
	saveFileFormat = id;
	int old_w = plr->width(), old_h = plr->height();
	int w = d.picWidth(), h = d.picHeight();
	if( w != old_w ||  h != old_h) {
	    plr->resize( w, h);
	    save();
	    plr->resize( old_w, old_h);
	} else
	    save();
    }
}


#include "plotQt.moc.cpp"
#endif // __FANCY_WIN__



//
// Implementation of the two architecture-dependent functions
// (from rect.h) requested by pari's plotting routines
//


void
rectdraw0(long *w, long *x, long *y, long lw)
{
    if (fork()) return;  // parent process returns

    pari_close();
    PARI_get_plot(1);

    // launch Qt window
    int argc = 1; char *argv[] = { "gp", "-qws"}; // set argc = 2 for cross
                                                  // development using qvfb
#ifdef __QPE__
    QPEApplication
#else
    QApplication
#endif
	a( argc, argv);
#ifdef __FANCY_WIN__
    PlotWindow *win = new PlotWindow(w, x, y, lw);
#else
    Plotter *win = new Plotter( w, x, y, lw);
#endif
#ifdef __QPE__
    a.showMainWidget( win);
#else
    a.setMainWidget( win);
    win->show();
#endif
    a.exec();
    exit( 0);
}

void
PARI_get_plot(long f)
/* This function initialises the structure rect.h: pari_plot */
{
    (void)f;
    if (pari_plot.init) return;      // pari_plot is already set
#ifdef __QPE__
    pari_plot.width   = 240;         // width and
    pari_plot.height  = 320;         //  height of plot window
#else
    pari_plot.width   = 400;         // width and
    pari_plot.height  = 300;         //  height of plot window
#endif
    pari_plot.hunit   = 3;           // 
    pari_plot.vunit   = 3;           //
    pari_plot.fwidth  = 6;           // font width
    pari_plot.fheight = 9;           //   and height
    pari_plot.init    = 1;           // flag: pari_plot is set now!
}
