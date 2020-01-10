#include "gifwriter.h"
#include <QImage>
#include <vector>
#include <QFile>
#include <QPainter>
#include <QLinearGradient>
#include "../gif.h"

int main()
{
  std::vector<QImage> imgs;
  for (int i = 0; i < 2; ++i) {
    QImage img(660, 660, MyLib::GifWriter::image_format);
    QPainter painter(&img);

    {
      QLinearGradient grad;
      grad.setCoordinateMode(QGradient::ObjectMode);
      grad.setStart(QPointF(0.0, 0.0));
      grad.setFinalStop(QPointF(1.0, 1.0));
      grad.setColorAt(0.0,  QColor(Qt::red));
      grad.setColorAt(0.25, QColor(Qt::blue));
      grad.setColorAt(0.5,  QColor(Qt::yellow));
      grad.setColorAt(0.75, QColor(Qt::green));
      grad.setColorAt(1.0,  QColor(Qt::white));
      painter.fillRect(img.rect(), grad);
    }

    {
      QLinearGradient grad;
      grad.setCoordinateMode(QGradient::ObjectMode);
      grad.setStart(QPointF(1.0, 0.0));
      grad.setFinalStop(QPointF(0.0, 1.0));
      grad.setColorAt(i == 0 ? 0.0 : 1.0, Qt::white);
      grad.setColorAt(i == 0 ? 1.0 : 0.0, Qt::black);
      painter.setCompositionMode(QPainter::CompositionMode_Multiply);
      painter.fillRect(img.rect(), grad);
    }

    imgs.push_back(img);
  }

  int delay = 200;
  bool dither = true;

  // use the new cpp-implementation
  QFile file("/tmp/test-cpp.gif");
  if (file.open(QIODevice::WriteOnly)) {
    const auto data = MyLib::GifWriter::encode(imgs, delay, dither);
    file.write(data);
  }

  // use the old base implementation to verify the results.
  GifWriter writer;
  GifBegin(&writer, "/tmp/test-old.gif", imgs.front().width(), imgs.front().height(), delay, 8, dither);
  for (const QImage& frame : imgs) {
    GifWriteFrame(&writer, frame.bits(), frame.width(), frame.height(), delay, 8, dither);
  }
  GifEnd(&writer);
}
