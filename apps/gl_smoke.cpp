#include <QApplication>
#include <QElapsedTimer>
#include <QOpenGLBuffer>
#include <QOpenGLFunctions>
#include <QOpenGLShaderProgram>
#include <QOpenGLVertexArrayObject>
#include <QOpenGLWidget>
#include <QSurfaceFormat>
#include <QTimer>
#include <QVector4D>
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <vector>

// Draws 2 vertices per device-pixel column (a min/max style envelope), every frame, and reports FPS.
class GlProbe : public QOpenGLWidget, protected QOpenGLFunctions {
public:
    GlProbe() : vbo_(QOpenGLBuffer::VertexBuffer) {
        setWindowTitle("RadarSim GL probe");
        resize(1000, 500);
    }

protected:
    void initializeGL() override {
        initializeOpenGLFunctions();
        auto s = [this](GLenum e) { return reinterpret_cast<const char*>(glGetString(e)); };
        const QSurfaceFormat f = format();
        std::printf("GL_VENDOR   : %s\n", s(GL_VENDOR));
        std::printf("GL_RENDERER : %s\n", s(GL_RENDERER));
        std::printf("GL_VERSION  : %s\n", s(GL_VERSION));
        std::printf("GLSL        : %s\n", s(GL_SHADING_LANGUAGE_VERSION));
        std::printf("Context     : %d.%d %s profile, swap interval %d\n", f.majorVersion(), f.minorVersion(),
                    f.profile() == QSurfaceFormat::CoreProfile ? "core" : "non-core", f.swapInterval());
        const char* r = s(GL_RENDERER);
        if (r && (std::strstr(r, "llvmpipe") || std::strstr(r, "softpipe") || std::strstr(r, "SWR")))
            std::puts("WARNING     : software rendering (no GPU). Expect low FPS at large sizes.");

        const bool ok = prog_.addShaderFromSourceCode(QOpenGLShader::Vertex,
                            "#version 330 core\nlayout(location=0) in vec2 pos;\nvoid main(){ gl_Position = vec4(pos, 0.0, 1.0); }") &&
                        prog_.addShaderFromSourceCode(QOpenGLShader::Fragment,
                            "#version 330 core\nuniform vec4 color;\nout vec4 frag;\nvoid main(){ frag = color; }") &&
                        prog_.link();
        std::printf("Shaders     : %s\n", ok ? "compiled and linked" : "FAILED");
        if (!ok) std::printf("%s\n", prog_.log().toUtf8().constData());
        vao_.create();
        vbo_.create();

        auto* t = new QTimer(this);
        QObject::connect(t, &QTimer::timeout, this, [this] { update(); });
        t->start(0);
        clock_.start();
    }

    void paintGL() override {
        glClearColor(0.07f, 0.08f, 0.10f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT);
        const int cols = std::max(2, static_cast<int>(width() * devicePixelRatioF()));
        verts_.resize(static_cast<std::size_t>(cols) * 4);
        phase_ += 0.03;
        for (int c = 0; c < cols; ++c) {
            const double u = static_cast<double>(c) / cols;
            const float x = static_cast<float>(-1.0 + 2.0 * (c + 0.5) / cols);
            const double mid = 0.5 * std::sin(6.2831853 * 3.0 * u + phase_);
            const double half = 0.15 * (1.2 + std::sin(6.2831853 * 11.0 * u - phase_));
            verts_[4 * c + 0] = x; verts_[4 * c + 1] = static_cast<float>(mid - half);
            verts_[4 * c + 2] = x; verts_[4 * c + 3] = static_cast<float>(mid + half);
        }
        vao_.bind();
        prog_.bind();
        prog_.setUniformValue("color", QVector4D(0.30f, 0.80f, 1.0f, 1.0f));
        vbo_.bind();
        vbo_.allocate(verts_.data(), static_cast<int>(verts_.size() * sizeof(float)));
        prog_.enableAttributeArray(0);
        prog_.setAttributeBuffer(0, GL_FLOAT, 0, 2);
        glDrawArrays(GL_LINES, 0, static_cast<GLsizei>(verts_.size() / 2));
        vbo_.release();
        prog_.release();
        vao_.release();

        ++frames_;
        if (clock_.elapsed() >= 1000) {
            const double fps = frames_ * 1000.0 / static_cast<double>(clock_.elapsed());
            std::printf("FPS %.1f  (%d columns, %zu vertices per frame)\n", fps, cols, verts_.size() / 2);
            setWindowTitle(QString("RadarSim GL probe  -  %1 FPS").arg(fps, 0, 'f', 1));
            frames_ = 0;
            clock_.restart();
        }
    }

private:
    QOpenGLShaderProgram prog_;
    QOpenGLBuffer vbo_;
    QOpenGLVertexArrayObject vao_;
    std::vector<float> verts_;
    QElapsedTimer clock_;
    double phase_ = 0.0;
    int frames_ = 0;
};

int main(int argc, char** argv) {
    QSurfaceFormat fmt;
    fmt.setVersion(3, 3);
    fmt.setProfile(QSurfaceFormat::CoreProfile);
    fmt.setSwapInterval(1);
    QSurfaceFormat::setDefaultFormat(fmt);

    QApplication app(argc, argv);
    GlProbe w;
    w.show();
    if (argc > 1) QTimer::singleShot(static_cast<int>(std::atof(argv[1]) * 1000.0), &app, &QCoreApplication::quit);
    return app.exec();
}
