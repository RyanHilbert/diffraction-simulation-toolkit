#ifndef THREAD_H
#define THREAD_H
#include<QThread>
class Thread:public QThread{
	Q_OBJECT
public:
	using QThread::QThread;
	void run()override;
public slots:
	void begin();//emits the starting() signal and starts the thread
	void stop();
signals:
	void starting();//emitted from the calling thread, not the new thread, unlike started()
	void error(QString str);
	void warning(QString str);
	void output(QString str);
	void maximum(int max);
	void progress(int prg);
	void completed1st(int garbage=0);
};
extern Thread calculation;
#endif