// CMakeProject1.cpp: определяет точку входа для приложения.
//

#include "pendulum_sim.h"

#include <Mahi/Gui.hpp>
#include <Mahi/Util.hpp>
#include <string>

struct ScrollingBuffer {
	int MaxSize;
	int Offset;
	ImVector<ImVec2> Data;
	ScrollingBuffer(int max_size = 2000) {
		MaxSize = max_size;
		Offset = 0;
		Data.reserve(MaxSize);
	}
	void AddPoint(float x, float y) {
		if (Data.size() < MaxSize)
			Data.push_back(ImVec2(x, y));
		else {
			Data[Offset] = ImVec2(x, y);
			Offset = (Offset + 1) % MaxSize;
		}
	}
	void Erase() {
		if (Data.size() > 0) {
			Data.shrink(0);
			Offset = 0;
		}
	}
};

class PendulumApp : public mahi::gui::Application {
private:
	mahi::gui::Vec2 worldOrigin = { 0.0f, 0.0f }; //world origin in screen space
	float scale = 100; //pixels in one meter
	float infoWindowSize = 450;

	double simulationTime = 0;

	//PHYSICS VARIABLES
	double g = 9.81;
    double l = 1;
	double lambda = 0;
	double fixedDt = 0.01666;

	//INITIAL VALUES
	double initial_fi = 1.57;

	double acc = 0;
	double vel = 0;
	double fi = initial_fi;

	bool started = false;
	bool simulating = false;
	bool useVerlet = true;

	//CALCULATED PERIOD
	double huygens = 0;
	double CEI = 0;
	double measured = 0;

	//FOR PERIOD MEASUREMENT
	double temp = 0;
	double t0 = 0;
	int Nt = 0;

	//LOGGING
	ScrollingBuffer fi_data;
	ScrollingBuffer vel_data;
	ScrollingBuffer acc_data;
	ScrollingBuffer phase_data;

public:
	PendulumApp() : Application(1200, 600, "pendulum_sim") {
		glfwSwapInterval(1);
		on_window_resized.connect(this, &PendulumApp::window_resize_handler);
		
		worldOrigin.x = (get_window_size().x - infoWindowSize) / 2.f;
		worldOrigin.y = get_window_size().y / 2.f;
	}

	void Euiler() {
		acc = ((-(g / l)) * std::sin(fi)) - 2 * lambda * vel;
		vel += acc * fixedDt;
		fi = fi + vel * fixedDt;
	}

	void Verlet() {
		double vel_pr, acc_pr;
		acc = ((-(g / l)) * std::sin(fi)) - 2 * lambda * vel;
		fi = fi + vel * fixedDt + 0.5 * acc * fixedDt * fixedDt;
		vel_pr = vel + acc * fixedDt;
		acc_pr = ((-(g / l)) * std::sin(fi)) - 2 * lambda * vel_pr;
		vel = vel + 0.5 * (acc + acc_pr) * fixedDt;
	}

	void update() override {
		//for(int i = 0; i < 4; i++)
		if (simulating) {
			double fi_last = fi;
			if (useVerlet) Verlet();
			else Euiler();

			if (fi >= 0 && fi_last < 0) {
				temp = fixedDt * fi_last / (fi_last - fi);
				if (Nt == 0) {
					t0 = simulationTime + temp;
				}
				else {
					measured = (simulationTime + temp - t0) / Nt;
				}
				Nt++;
			}

			simulationTime += fixedDt;


			fi_data.AddPoint(simulationTime, fi);
			vel_data.AddPoint(simulationTime, vel);
			acc_data.AddPoint(simulationTime, acc);
			phase_data.AddPoint(fi, vel);
		}

		huygens = 2.0 * 3.1415926 * std::sqrt(l / g);
		CEI = 4.0 * std::sqrt(l / g) * cei1(std::sin(initial_fi / 2.0));
		
		if (!started) fi = initial_fi;

		drawPendulum();

		

		//UI
		ImVec2 pos = ImGui::GetMainViewport()->Pos;
		ImVec2 size = ImGui::GetMainViewport()->Size;

		ImGui::SetNextWindowPos({ pos.x + size.x - infoWindowSize, pos.y });
		ImGui::SetNextWindowSize({ infoWindowSize, size.y});
		ImGui::Begin("Info window", NULL, ImGuiWindowFlags_NoCollapse | ImGuiWindowFlags_NoMove | ImGuiWindowFlags_NoResize | ImGuiWindowFlags_NoDocking | ImGuiWindowFlags_NoTitleBar);
		ImGui::Text("FPS: %.3f", ImGui::GetIO().Framerate);
		if (ImGui::Button("Play", ImVec2(50, 0))) {
			simulating = true;
			started = true;
		}
		ImGui::SameLine();
		if (ImGui::Button("Pause", ImVec2(50, 0))) {
			simulating = false;
		}
		ImGui::SameLine();
		if (ImGui::Button("Reset", ImVec2(50, 0))) {
			started = false;
			simulating = false;
			fi = initial_fi;
			vel = 0;
			acc = 0;
			simulationTime = 0;
			temp = 0;
			Nt = 0;
			t0 = 0;
			fi_data.Erase();
			vel_data.Erase();
			acc_data.Erase();
			phase_data.Erase();
		}

		ImGui::Checkbox("Use Verlet", &useVerlet);
		ImGui::Text("Variables:");
		
		ImGui::Indent();
		ImGui::SliderDouble("Initial fi", &initial_fi, 0, 6.2831853, "%.3f radians");
		ImGui::SliderDouble("dt", &fixedDt, 0.001, 0.5, "%.3f s");
		ImGui::SliderDouble("G", &g, 0, 30, "%.3f m/s^2");
		ImGui::SliderDouble("Lambda", &lambda, 0, 3);
		ImGui::SliderDouble("l", &l , 0.1, 5);
		ImGui::Unindent();

		ImGui::Spacing();
		ImGui::Text("Calculated period:");

		std::string sHuygens = "Huygens: " + std::to_string(huygens) + " s.";
		ImGui::BulletText(sHuygens.c_str());

		std::string sCEI = "CEI: " + ((lambda < 0.0001) ? std::to_string(CEI) + " s." : "damped oscillation.");
		ImGui::BulletText(sCEI.c_str());

		std::string sMeasured = "Measured: " + std::to_string(measured) + " s.";
		ImGui::BulletText(sMeasured.c_str());

		ImGui::Spacing();

		static float history = 10.f;
		ImGui::SliderFloat("Plot history", &history, 1, 30, "%.1f s");

		ImPlot::SetNextPlotLimitsX(simulationTime - history, simulationTime, ImGuiCond_Always);
		ImPlot::FitNextPlotAxes(false, true, false, false);
		if (ImPlot::BeginPlot("FI plot", "time (s)", "angle (rad)")) {
			if(!fi_data.Data.empty())
				ImPlot::PlotLine("fi", &fi_data.Data[0].x, &fi_data.Data[0].y, fi_data.Data.size(), fi_data.Offset, 2 * sizeof(float));
			if (!vel_data.Data.empty())
				ImPlot::PlotLine("vel", &vel_data.Data[0].x, &vel_data.Data[0].y, vel_data.Data.size(), vel_data.Offset, 2 * sizeof(float));
			if (!acc_data.Data.empty())
				ImPlot::PlotLine("acc", &acc_data.Data[0].x, &acc_data.Data[0].y, acc_data.Data.size(), acc_data.Offset, 2 * sizeof(float));
			ImPlot::EndPlot();
		}

		static bool autoScale = true;
		ImGui::Checkbox("Auto fit phase portrait", &autoScale);
		if(autoScale) ImPlot::FitNextPlotAxes(true, true, false, false);
		if (ImPlot::BeginPlot("Phase portrait", "angle (rad)", "angular velocity (rad/s)")) {
			if (!fi_data.Data.empty())
				ImPlot::PlotLine("phase", &phase_data.Data[0].x, &phase_data.Data[0].y, phase_data.Data.size(), phase_data.Offset, 2 * sizeof(float));
			ImPlot::EndPlot();
		}

		ImGui::End();

		

	}

private:

	double cei1(double k) {
		double a, b, t;
		t = 1 - k * k;
		a = (((0.01451196212 * t + 0.03742563713) * t
			+ 0.03590092383) * t + 0.09666344259) * t + 1.38629436112;
		b = (((0.00441787012 * t + 0.03328355346) * t + 0.06880248576)
			* t + 0.12498593597) * t + 0.5;
		return a - b * std::log(t);
	}

	void drawPendulum() {
		mahi::gui::Vec2 pos = { (float)(std::sin(fi) * l), (float)(std::cos(fi) * l) };
		mahi::gui::Vec2 screenPos = worldToScreen(pos);

		nvgBeginPath(m_vg);
		nvgStrokeColor(m_vg, nvgRGBA(100, 100, 100, 255));
		nvgStrokeWidth(m_vg, 1);

		nvgEllipse(m_vg, screenPos.x, screenPos.y, 5.f, 5.f);

		nvgMoveTo(m_vg, worldOrigin.x, worldOrigin.y);
		nvgLineTo(m_vg, screenPos.x, screenPos.y);

		nvgStroke(m_vg);
	}

	void drawGrid() {
		nvgBeginPath(m_vg);
		nvgStrokeColor(m_vg, nvgRGBA(100,100,100,255));
		nvgStrokeWidth(m_vg, 1);

		nvgEllipse(m_vg, worldOrigin.x, worldOrigin.y, 5.f, 5.f);

		nvgMoveTo(m_vg, 0, 10);
		nvgLineTo(m_vg, get_window_size().x, 10);
		nvgMoveTo(m_vg, 0, 20);
		nvgLineTo(m_vg, get_window_size().x, 20);
		nvgStroke(m_vg);
	}

	mahi::gui::Vec2 screenToWorld(mahi::gui::Vec2 pos) {
		return (pos - worldOrigin) / scale;
	}

	mahi::gui::Vec2 worldToScreen(mahi::gui::Vec2 pos) {
		return worldOrigin + pos * scale;
	}

	void window_resize_handler(int width, int height) {
		worldOrigin.x = (get_window_size().x - infoWindowSize) / 2.f;
		worldOrigin.y = get_window_size().y / 2.f;
	}

};

int main() {
	PendulumApp app;
	app.run();
	return 0;
}
