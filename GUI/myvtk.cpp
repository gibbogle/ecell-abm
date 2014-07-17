// myvtk.cpp

#include <vtkCamera.h>
#include <vtkInteractorStyleTrackballCamera.h>
#include <vtkObjectFactory.h>

#ifdef _WIN32
#include "windows.h"
#endif
#include "myvtk.h"
#include "log.h"
#include "transfer.h"

LOG_USE();

// Define interaction style
class MouseInteractorStyle4 : public vtkInteractorStyleTrackballCamera
{
  public:
	static MouseInteractorStyle4* New();
	vtkTypeMacro(MouseInteractorStyle4, vtkInteractorStyleTrackballCamera);

	virtual void OnLeftButtonDown()
	{
	  leftb = true;
	  // Forward events
	  vtkInteractorStyleTrackballCamera::OnLeftButtonDown();
	}

	virtual void OnMiddleButtonDown()
	{
//	  std::cout << "Pressed middle mouse button." << std::endl;
	  // Forward events
	  vtkInteractorStyleTrackballCamera::OnMiddleButtonDown();
	}

	virtual void OnRightButtonDown()
	{
//	  std::cout << "Pressed right mouse button." << std::endl;
	  // Forward events
	  vtkInteractorStyleTrackballCamera::OnRightButtonDown();
	}

	virtual void OnLeftButtonUp()
	{
//	  std::cout << "Released left mouse button." << std::endl;
//	  LOG_QMSG("Released left mouse button.");
	  leftb = false;
	  // Forward events
	  vtkInteractorStyleTrackballCamera::OnLeftButtonUp();
	}

	virtual void OnMiddleButtonUp()
	{
//	  std::cout << "Released middle mouse button." << std::endl;
	  // Forward events
	  vtkInteractorStyleTrackballCamera::OnMiddleButtonUp();
	}

	virtual void OnRightButtonUp()
	{
//	  std::cout << "Released right mouse button." << std::endl;
	  // Forward events
	  vtkInteractorStyleTrackballCamera::OnRightButtonUp();
	}

};

vtkStandardNewMacro(MouseInteractorStyle4);

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
MyVTK::MyVTK(QWidget *page, QWidget *key_page)
{
	zoomlevel = 0.7;
	double backgroundColor[] = {0.0,0.0,0.0};


	Pi = 4*atan(1.0);
    page_VTK = page;
	leftb = false;
    key_canvas(key_page);
    qvtkWidget = new QVTKWidget(page,QFlag(0));
	LOG_MSG("Created a new QVTKWidget");
	QVBoxLayout *layout = new QVBoxLayout;
    layout->addWidget(qvtkWidget);

	// Associate the layout with page_VTK
    page_VTK->setLayout(layout);

	// Create a renderer, and add it to qvtkWidget's render window.
	// The renderer renders into the render window. 
    ren = vtkSmartPointer<vtkRenderer>::New();
    renWin = qvtkWidget->GetRenderWindow();
    renWin->AddRenderer(ren);
	ren->SetBackground(backgroundColor);
	ren->ResetCamera();
	iren = qvtkWidget->GetInteractor();

	vtkSmartPointer<MouseInteractorStyle4> style = vtkSmartPointer<MouseInteractorStyle4>::New();
	iren->SetInteractorStyle( style );

	iren->Initialize();

	// Create mappers
    createMappers();
    Eactor = vtkSmartPointer<vtkActor>::New();
    Eactor->SetMapper(Emapper);

// Create image filter for save Snapshot()
//	w2img = vtkWindowToImageFilter::New();
//	pngwriter = vtkSmartPointer<vtkPNGWriter>::New();
//	jpgwriter = vtkSmartPointer<vtkJPEGWriter>::New();

	first_VTK = true;
	DCmotion = false;
    DCfade = false;
	playing = false;
	paused = false;
    opacity = 1.0;
    display_celltype[1] = true;
    display_celltype[2] = true;
    TCpos_list.clear();
    ren->GetActiveCamera()->Zoom(zoomlevel);		// try zooming OUT
}

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
MyVTK::~MyVTK()
{
}

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
void MyVTK::key_canvas(QWidget *key_page)
{
//    QGraphicsScene* scene = new QGraphicsScene(QRect(0, 0, 130, 280));
    QGraphicsScene* scene = new QGraphicsScene(QRect(0, 0, 130, 310));
    QBrush brush;
    QGraphicsTextItem *text;

//    brush.setColor(QColor(255,128,77));
    brush.setColor(QColor(150,100,0));      // dark brown
    brush.setStyle(Qt::SolidPattern);
	scene->addEllipse(10,10,20,20,Qt::NoPen, brush);
	text = scene->addText("FDC");
    text->setPos(35, 10);

//    brush.setColor(QColor(255,77,128));
    brush.setColor(QColor(200,60,100));     // dark red
    brush.setStyle(Qt::SolidPattern);
    scene->addEllipse(10,40,20,20,Qt::NoPen, brush);
    text = scene->addText("MRC");
    text->setPos(35, 40);

    brush.setColor(QColor(30,20,255));      // dark blue
    scene->addEllipse(10,70,20,20,Qt::NoPen, brush);
	text = scene->addText("Naive B cell");
    text->setPos(35, 70);

    brush.setColor(QColor(0,200,255));      // light blue
    scene->addEllipse(10,100,20,20,Qt::NoPen, brush);
	text = scene->addText("CCR7 UP");
    text->setPos(35, 100);

    brush.setColor(QColor(50,255,150));     // light green
    scene->addEllipse(10,130,20,20,Qt::NoPen, brush);
	text = scene->addText("EBI2 UP");
    text->setPos(35, 130);

//    brush.setColor(QColor(255,255,0));      // yellow
    brush.setColor(Qt::yellow );      // yellow
    scene->addEllipse(10,160,20,20,Qt::NoPen, brush);
	text = scene->addText("BCL6 HI");
    text->setPos(35, 160);

    brush.setColor(QColor(0,150,0));        // dark green
    scene->addEllipse(10,190,20,20,Qt::NoPen, brush);
	text = scene->addText("BCL6 LO");
    text->setPos(35, 190);

    brush.setColor(QColor(128,128,128));    // grey
    scene->addEllipse(10,220,20,20,Qt::NoPen, brush);
	text = scene->addText("Max divisions");
    text->setPos(35, 220);

    brush.setColor(QColor(255,0,0));        // red
    scene->addEllipse(10,250,20,20,Qt::NoPen, brush);
	text = scene->addText("Plasma cell");
    text->setPos(35, 250);

    brush.setColor(QColor(255,130,0));      // orange
    scene->addEllipse(10,280,20,20,Qt::NoPen, brush);
    text = scene->addText("CD4 T cell");
    text->setPos(35, 280);

	QGraphicsView* view = new QGraphicsView(key_page);
    view->setScene(scene);
//    view->setGeometry(QRect(0, 0, 150, 300));
    view->setGeometry(QRect(0, 0, 150, 330));
    view->show();
}

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
void MyVTK::createMappers()
{

// Ellipsoid mapper
    Epoints = vtkSmartPointer<vtkPoints>::New();
    EpolyData = vtkSmartPointer<vtkPolyData>::New();
    EpolyData->SetPoints(Epoints);
    Etensors = vtkSmartPointer<vtkDoubleArray>::New();
    Etensors->SetNumberOfComponents(9);
    EpolyData->GetPointData()->SetTensors(Etensors);

    EsphereSource = vtkSmartPointer<vtkSphereSource>::New();
    EsphereSource->Update();
    EsphereSource->SetThetaResolution(16);
    EsphereSource->SetPhiResolution(16);
    EsphereSource->SetRadius(1.0);

    EtensorGlyph = vtkSmartPointer<vtkTensorGlyph>::New();
//    EtensorGlyph->SetSourceConnection(EsphereSource->GetOutputPort());   // Both methods work
    EtensorGlyph->SetSource(EsphereSource->GetOutput());
    EtensorGlyph->SetInput(EpolyData);
    EtensorGlyph->ColorGlyphsOff();
    EtensorGlyph->ThreeGlyphsOff();
    EtensorGlyph->ExtractEigenvaluesOff();
//    EtensorGlyph->Update();

    Emapper = vtkSmartPointer<vtkPolyDataMapper>::New();
//    Emapper->SetInputConnection(EtensorGlyph->GetOutputPort());   // Both methods work
    Emapper->SetInput(EtensorGlyph->GetOutput());
}

//-----------------------------------------------------------------------------------------
// The cell info is fetched from the DLL by ExecThread::snapshot().
// The info is transmitted in the double array ecell_list[],
// each cell having 3 position values and 9 tensor values
//-----------------------------------------------------------------------------------------
void MyVTK::get_cell_positions()
{
    LOG_QMSG("get_cell_positions");
    for (int i=0; i<necell_list; i++) {
        int j = NINFO*i;
        Epoints->InsertPoint(i,
                  ecell_list[j],  ecell_list[j+1], ecell_list[j+2]);
        Etensors->InsertTuple9(i,
                  ecell_list[j+3],ecell_list[j+4], ecell_list[j+5],
                  ecell_list[j+6],ecell_list[j+7], ecell_list[j+8],
                  ecell_list[j+9],ecell_list[j+10],ecell_list[j+11]);
	}
}

//---------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------
void MyVTK::init()
{
    Epoints->Reset();
    Etensors->Reset();
//    EpolyData->Reset();
//    EpolyData->SetPoints(Epoints);
//    EpolyData->GetPointData()->SetTensors(Etensors);
//    EtensorGlyph->SetInput(EpolyData);
//    Emapper = vtkSmartPointer<vtkPolyDataMapper>::New();
//    Emapper->SetInput(EtensorGlyph->GetOutput());
//    Eactor->SetMapper(Emapper);
    sprintf(msg,"init: EpolyData->GetNumberOfPoints(): %d",EpolyData->GetNumberOfPoints());
    LOG_MSG(msg);
}

//---------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------
void MyVTK::cleanup()
{
	LOG_MSG("VTK cleanup");
	first_VTK = true;	
}

//---------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------
void MyVTK::renderCells(bool redo, bool zzz)
{
    LOG_QMSG("renderCells");
    process_Ecells();
	if (first_VTK) {
		LOG_MSG("Initializing the renderer");
		ren->ResetCamera();
	}
	iren->Render();
    first_VTK = false;
//    if (record) {
//        recorder();
//    }
}

//---------------------------------------------------------------------------------------------
// Interprets an int as rgb
// USE_CELLTYPE_COLOUR, the cell type is passed in cp.state, and the colours are those
// that were chosen in the GUI.
// Otherwise cp.state is packed (r,g,b)
//---------------------------------------------------------------------------------------------
void MyVTK::unpack(int x, double *rr, double *gg, double *bb)
{
	int z, r, g, b;

    if (USE_CELLTYPE_COLOUR) {

    } else {
        z = x;
        r = z>>16;
        z = r;
        z = z<<16;

        x = x - z;

        z = x;
        g = z>>8;
        z = g;
        z = z<<8;

        b = x - z;
    }
    *rr = r/255.;
    *gg = g/255.;
    *bb = b/255.;
}

//---------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------
void MyVTK::process_Ecells()
{
    EtensorGlyph->Update();
    EtensorGlyph->RemoveAllInputs();
    EtensorGlyph->SetInput(EpolyData);
    int na = ren->GetActors()->GetNumberOfItems();
    if (na == 0) {
        LOG_QMSG("AddActor");
        ren->AddActor(Eactor);
    }
}

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
bool MyVTK::startPlayer(QString posfile, QTimer *theTimer, bool save)
{
	save_image = save;
	LOG_QMSG(posfile);
	timer = theTimer;
	playerData = new QFile(posfile);
	if (!playerData->open(QFile::ReadOnly)) {
		LOG_MSG("Open failure on VTK file");
		return false;
	}
	playerStream = new QTextStream(playerData);
	if (!first_VTK) {
		cleanup();
	}
	playing = true;
	paused = false;

	if (save_image) {
        w2i = vtkWindowToImageFilter::New();
        w2i->SetInput(renWin);	//the render window
//		writer = vtkSmartPointer<vtkPNGWriter>::New();
        jpgwriter = vtkSmartPointer<vtkJPEGWriter>::New();
        jpgwriter->SetInputConnection(w2i->GetOutputPort());
		framenum = 0;
		LOG_MSG("set up writer");
	}
	LOG_MSG("playing");
	return true;
}

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
bool MyVTK::nextFrame()
{
	LOG_MSG("VTK: nextFrame");
	if (!playing)
		return false;
	if (paused)
		return true;
	if (playerStream->atEnd()) {
		LOG_MSG("nextFrame: no more data");
		stop();
		return false;
	}
//    TCpos_list.clear();
//	DCpos_list.clear();
//	bondpos_list.clear();
	int k = 0;
	QString line;
	do {
		line = playerStream->readLine();
		if (line.length() > 0) {
			k++;
			QStringList s = line.split(" ",QString::SkipEmptyParts);
			if (s[0].compare("T") == 0) {
				CELL_POS cp;
				cp.tag = s[1].toInt();
				cp.x = s[2].toInt();
				cp.y = s[3].toInt();
				cp.z = s[4].toInt();
				cp.diameter = s[5].toDouble();
				cp.state = s[6].toDouble();
                TCpos_list.append(cp);
			} else if (s[0].compare("E") == 0) {
				break;
			}
		}
	} while (true);

	bool redo = false;
	if (first_VTK) {
		redo = true;
	}
    renderCells(redo,false);
	char numstr[5];
	sprintf(numstr,"%04d",framenum);
	if (save_image) {
        w2i->Modified();	//important
        jpgwriter->SetFileName((casename + numstr + ".jpg").toStdString().c_str());
        jpgwriter->Write();
	}
	framenum++;
	return true;
}

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
void MyVTK::pause()
{
    paused = true;
}

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
void MyVTK::playon()
{
    paused = false;
}

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
void MyVTK::stop()
{
    if (save_image) {
        if (jpgwriter) jpgwriter->Delete();
        if (pngwriter) pngwriter->Delete();
        w2i->Delete();
    }
    delete playerStream;
    playerData->close();
    delete playerData;
    timer->stop();
    playing = false;
    paused = false;
}


//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
void MyVTK::saveSnapshot(QString fileName, QString imgType)
{
    vtkSmartPointer<vtkWindowToImageFilter> w2img = vtkWindowToImageFilter::New();

	w2img->SetInput(renWin);
	if (imgType.compare("png") == 0) {
		vtkSmartPointer<vtkPNGWriter> pngwriter = vtkPNGWriter::New();
		pngwriter->SetInputConnection(w2img->GetOutputPort()); 
		w2img->Modified();
		pngwriter->SetFileName((fileName.toStdString()).c_str()); 
		pngwriter->Write();
//		pngwriter->Delete();	// Note: using vtkSmartPointer, delete is not necessary.
	} else if (imgType.compare("jpg") == 0) {
		vtkJPEGWriter *jpgwriter = vtkJPEGWriter::New();
		jpgwriter->SetInputConnection(w2img->GetOutputPort()); 
		w2img->Modified();
		jpgwriter->SetFileName((fileName.toStdString()).c_str()); 
		jpgwriter->Write();
//		jpgwriter->Delete();
	} else if (imgType.compare("tif") == 0) {
		vtkTIFFWriter *tifwriter = vtkTIFFWriter::New();
		tifwriter->SetInputConnection(w2img->GetOutputPort()); 
		w2img->Modified();
		tifwriter->SetFileName((fileName.toStdString()).c_str()); 
		tifwriter->Write();
//		tifwriter->Delete();
	} else if (imgType.compare("bmp") == 0) {
		vtkBMPWriter *bmpwriter = vtkBMPWriter::New();
		bmpwriter->SetInputConnection(w2img->GetOutputPort()); 
		w2img->Modified();
		bmpwriter->SetFileName((fileName.toStdString()).c_str()); 
		bmpwriter->Write();
	}
}


//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
void MyVTK::set_celltype_colour(COLOUR_TYPE *colour, QString str)
{
    if (str.compare("red") == 0) {
        colour->r = 1.0;
        colour->g = 0.0;
        colour->b = 0.0;
    } else if (str.compare("orange") == 0) {
        colour->r = 0.8;
        colour->g = 0.5;
        colour->b = 0.0;
    } else if (str.compare("yellow") == 0) {
        colour->r = 1.0;
        colour->g = 1.0;
        colour->b = 0.0;
    } else if (str.compare("green") == 0) {
        colour->r = 0.0;
        colour->g = 1.0;
        colour->b = 0.0;
    } else if (str.compare("blue") == 0) {
        colour->r = 0.0;
        colour->g = 0.0;
        colour->b = 1.0;
    } else if (str.compare("purple") == 0) {
        colour->r = 1.0;
        colour->g = 0.0;
        colour->b = 1.0;
    } else if (str.compare("brown") == 0) {
        colour->r = 0.5;
        colour->g = 0.5;
        colour->b = 0.2;
    }

}

/*
//-----------------------------------------------------------------------------------------
// Uses: w2i, renWin, pngwriter,
//       record, record_basename, record_nframes, record_it, framenum
//-----------------------------------------------------------------------------------------
void MyVTK::startRecorder(QString basefile, int nframes)
{
    if (w2i == 0) {
        w2i = vtkWindowToImageFilter::New();
        w2i->SetInput(renWin);	//the render window
//        jpgwriter = vtkSmartPointer<vtkJPEGWriter>::New();
//        jpgwriter->SetInputConnection(w2i->GetOutputPort());
        pngwriter = vtkSmartPointer<vtkPNGWriter>::New();
        pngwriter->SetInputConnection(w2i->GetOutputPort());
}
    framenum = 0;
    LOG_MSG("set up writer");
    record = true;
    record_basename = basefile;
    record_nframes = nframes;
    record_it = 0;

    LOG_MSG("Started recording");
}

//-----------------------------------------------------------------------------------------
// Uses: w2i, videoOutput, pngwriter,
//       record, record_basename, record_nframes, record_it, framenum
//-----------------------------------------------------------------------------------------
void MyVTK::recorder()
{
    char numstr[5];
    char filename[512];

    sprintf(msg,"recorder: record_it: %d",record_it);
    LOG_MSG(msg);
    if (record_it > record_nframes) {
        record = false;
        stopRecorder();
        return;
    }
    vtkImageData *id = vtkImageData::New();
    id = w2i->GetOutput();
    w2i->Modified();	//important
    id->Update();
    int width = id->GetDimensions()[0];
    int height = id->GetDimensions()[1];
    if (width == 0) {
        LOG_QMSG("ERROR: recorder: vtkImageData dimension = 0");
        exit(1);
    }
    if (RECORD_VIDEO) {
        if (!videoOutput->isOpen()) {
            // Generate temporary filename
            tempFile = new QTemporaryFile("qt_temp.XXXXXX.avi");
            if (tempFile->open())
            {
               // Open media file and prepare for recording
               QString fileName = tempFile->fileName();
                bool recording = videoOutput->openMediaFile(width, height, fileName.toAscii().data());
                if (!recording) {
                    LOG_QMSG("ERROR: openMediaFile failed");
                    record = false;
                    return;
                }
            }
        }
        bool success = videoOutput->newVtkFrame(id);
        if (!success) {
            LOG_QMSG("ERROR: newVtkFrame failed");
            record = false;
            exit(1);
        }
        record_it++;
        return;
    }
//    strcpy(filename,record_basename);
//    strcat(filename,numstr);
//    strcat(filename,".jpg");
//    strcpy(filename ,(record_basename + numstr + ".jpg").toStdString().c_str());
//    jpgwriter->SetFileName(filename);
//    jpgwriter->Write();
    sprintf(numstr,"%05d",framenum);
    strcpy(filename ,(record_basename + numstr + ".png").toStdString().c_str());
    pngwriter->SetFileName(filename);
    pngwriter->Write();
    sprintf(msg,"recorder: it: %d frame: %d filename: %s  id dimensions: %d %d",record_it,framenum,filename,id->GetDimensions()[0],id->GetDimensions()[1]);
    LOG_MSG(msg);
    framenum++;
    record_it++;
}

//-----------------------------------------------------------------------------------------
// Uses: record
//-----------------------------------------------------------------------------------------
void MyVTK::stopRecorder()
{
    record = false;
    if (RECORD_VIDEO) {
        videoOutput->closeMediaFile();
//        QString fileName = QFileDialog::getSaveFileName(page_VTK,
//                                                        "Save File",    //tr("Save File"),
//                                                        QString(),
//                                                        "Videos (*.avi)");  //tr("Videos (*.avi)"));
        QString fileName = "testfile.avi";
        if (fileName.isNull() == false) {
           QFile::copy(tempFile->fileName(), fileName);
        }
        delete tempFile;
        tempFile = 0x0;
    }
//    jpgwriter->RemoveAllInputs();
//    jpgwriter->Delete();
//    w2i->RemoveAllInputs();
//    w2i->Delete();
    LOG_MSG("Stopped recording");
}
*/


/*
void MyVTK::recordingSlot()
{
   if (recording)
   {
      videoOutput->closeMediaFile();
      recordingButton->setIcon(QIcon(":/robotnavigator/resources/video-camera-png.png"));
      recordingButton->setToolTip("Record Video");
      recordingTimer->stop();
      recording = false;
      QString fileName = QFileDialog::getSaveFileName(this,
                                                      tr("Save File"),
                                                      QString(),
                                                      tr("Videos (*.avi)"));
      if (fileName.isNull() == false)
      {
         QFile::copy(tempFile->fileName(), fileName);
      }
      delete tempFile;
      tempFile = 0x0;
   }
   else
   {
      blinkCount = 0;
      // Generate temporary filename
      tempFile = new QTemporaryFile("qt_temp.XXXXXX.avi");
      if (tempFile->open())
      {
         // Open media file and prepare for recording
         QString fileName = tempFile->fileName();
         recording = videoOutput->openMediaFile(640, 480, fileName.toAscii().data());
      }
      // Change tool tip on video button and start blinking timer
      recordingButton->setToolTip("Stop Recording");
      recordingTimer->start(500);
   }
}
*/
