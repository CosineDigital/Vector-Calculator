#include "Walnut/Application.h"
#include "Walnut/EntryPoint.h"

#include "Walnut/Image.h"

#include <glm/glm.hpp>

#include <imgui.h>

// helper macros so each button has a unique identifier
#define COPY_BUTTON_STR(var) "Copy##"#var
#define PASTE_BUTTON_STR(var) "Paste##"#var

class Calculator {
public:
	virtual void Calculate() = 0;
	static constexpr float rad2deg = 180.f / 3.1415926535897932384626433832795f;
};

// data types to organize each tool
class PositionVector : public Calculator {
public:
	glm::vec3 P, Q, PQ;
	float PQmag;

	virtual void Calculate() {
		PQ = (Q - P);
		PQmag = std::sqrtf(PQ.x * PQ.x + PQ.y * PQ.y + PQ.z * PQ.z);
	}
};

class UnitVector : public Calculator {
public:
	glm::vec3 R, L, TL;
	float T, Rmag;

	virtual void Calculate() {

		Rmag = std::sqrtf(R.x * R.x + R.y * R.y + R.z * R.z);

		L = glm::normalize(R);

		TL = L * T;
	}
};

class Angles : public Calculator {
public:
	// input vector
	glm::vec3 R;
	float Rmag;
	// cosine angles
	float cos_x, cos_y, cos_z;
	// spherical coordinates, per ISO convention
	float s_radius, s_theta, s_azimuth;
	// cylindrical coordinates, per ISO convention
	float c_rho, c_azimuth, c_z;

	virtual void Calculate() {
		Rmag = sqrt(R.x * R.x + R.y * R.y + R.z * R.z);

		// cosine angles
		cos_x = std::acos(R.x / Rmag) * rad2deg;
		cos_y = std::acos(R.y / Rmag) * rad2deg;
		cos_z = std::acos(R.z / Rmag) * rad2deg;

		// spherical coordinates
		s_radius = Rmag;
		s_theta = std::acos(R.z / Rmag) * rad2deg;

		if (R.x > 0) {
			s_azimuth = std::atan(R.y / R.x) * rad2deg;
		}
		else if (R.x < 0 && R.y >= 0) {
			s_azimuth = std::atan(R.y / R.x) * rad2deg + 180.f;
		}
		else if (R.x < 0 && R.y < 0) {
			s_azimuth = std::atan(R.y / R.x) * rad2deg - 180.f;
		}
		else if (R.x == 0 && R.y > 0) {
			s_azimuth = 90;
		}
		else if (R.x == 0 && R.y < 0) {
			s_azimuth = -90;
		}
		else if (R.x == 0 && R.y == 0) {
			s_azimuth = std::numeric_limits<float>::infinity();
		}

		// cylindrical coordinates
		c_rho = s_radius * std::sin(s_theta / rad2deg);
		c_azimuth = s_azimuth;
		c_z = s_radius * std::cos(s_theta / rad2deg);
	}
};

class CrossProduct : public Calculator {
public:
	glm::vec3 r, F, M;
	float Mmag;

	virtual void Calculate() {
		M = glm::cross(r, F);

		Mmag = std::sqrtf(M.x * M.x + M.y * M.y + M.z * M.z);
	}
};

class DotProduct : public Calculator {
public:
	glm::vec3 P, Q;
	float R;

	virtual void Calculate() {
		R = glm::dot(P, Q);
	}
};

class Distance : public Calculator {
private:
	glm::vec3 d, n;
	float nmag, dmag, dprod;

public:
	glm::vec3 P, L0, L, PL0;
	float dist;

	virtual void Calculate() {
		PL0 = P - L0;

		dprod = glm::dot(PL0, L);

		if (dprod == 0) { // edge case, they are perpendicular
			dist = std::sqrt(PL0.x * PL0.x + PL0.y * PL0.y + PL0.z * PL0.z);
		}
		else {
			d = dprod * glm::normalize(L);
			n = glm::cross(PL0, d);

			nmag = std::sqrt(n.x * n.x + n.y * n.y + n.z * n.z);
			dmag = std::sqrt(d.x * d.x + d.y * d.y + d.z * d.z);

			dist = nmag / dmag;
		}
	}
};

class VCLayer : public Walnut::Layer
{
	static constexpr int disabledFlags = ImGuiInputTextFlags_ReadOnly | ImGuiSelectableFlags_Disabled;

	PositionVector pv;
	void DrawPositionVector() noexcept {
		ImGui::Begin("Position Vector");

		if (ImGui::InputFloat3("P##pv.P", &pv.P[0]))
			pv.Calculate();

		DrawCopyPaste(&pv, &pv.P, COPY_BUTTON_STR(pv.P), PASTE_BUTTON_STR(pv.P));

		if (ImGui::InputFloat3("Q##pv.Q", &pv.Q[0]))
			pv.Calculate();

		DrawCopyPaste(&pv, &pv.Q, COPY_BUTTON_STR(pv.Q), PASTE_BUTTON_STR(pv.Q));

		// Output

		ImGui::Text("");

		ImGui::InputFloat3("PQ##pv.PQ", &pv.PQ[0], "%.3f", disabledFlags);
		DrawCopy(&pv, &pv.PQ, COPY_BUTTON_STR(pv.PQ));

		ImGui::InputFloat("|PQ|##pv.PQmag", &pv.PQmag, 0, 0, "%.3f", disabledFlags);

		ImGui::End();
	}

	UnitVector uv;
	void DrawUnitVector() noexcept {
		ImGui::Begin("Unit Vector");

		if (ImGui::InputFloat3("R##uv.R", &uv.R[0]))
			uv.Calculate();

		DrawCopyPaste(&uv, &uv.R, COPY_BUTTON_STR(uv.R), PASTE_BUTTON_STR(uv.R));

		ImGui::InputFloat("|R|##uv.Rmag", &uv.Rmag, 0, 0, "%.3f", disabledFlags);

		// Output
		ImGui::Text("");

		if (ImGui::InputFloat3("Unit##uv.L", &uv.L[0], "%.3f", disabledFlags))
			uv.Calculate();

		DrawCopy(&uv, &uv.L, COPY_BUTTON_STR(uv.L));

		ImGui::Text("");

		if (ImGui::InputFloat("|T|", &uv.T))
			uv.Calculate();

		if (ImGui::InputFloat3("T along Unit", &uv.TL[0], "%.3f", disabledFlags))
			uv.Calculate();

		DrawCopy(&uv, &uv.TL, COPY_BUTTON_STR(uv.TL));

		ImGui::End();
	}

	Angles a;
	void DrawAngles() noexcept {

		ImGui::Begin("Angles");

		if (ImGui::InputFloat3("R##a.R", &a.R[0]))
			a.Calculate();

		DrawCopyPaste(&a, &a.R, COPY_BUTTON_STR(a.R), PASTE_BUTTON_STR(a.R));

		ImGui::Text("\nCosine Angles");

		ImGui::InputFloat("Theta x", &a.cos_x, 0, 0, "%.3f", disabledFlags);
		ImGui::InputFloat("Theta y", &a.cos_y, 0, 0, "%.3f", disabledFlags);
		ImGui::InputFloat("Theta z", &a.cos_z, 0, 0, "%.3f", disabledFlags);

		ImGui::Text("\nSpherical Coordinates");

		ImGui::InputFloat("Radius##a.s_radius", &a.s_radius, 0, 0, "%.3f", disabledFlags);
		ImGui::InputFloat("Theta##a.s_theta", &a.s_theta, 0, 0, "%.3f", disabledFlags);
		ImGui::InputFloat("Azimuth##a.s_azimuth", &a.s_azimuth, 0, 0, "%.3f", disabledFlags);

		ImGui::Text("\nCylindrical Coordinates");

		ImGui::InputFloat("Rho##a.c_rho", &a.c_rho, 0, 0, "%.3f", disabledFlags);
		ImGui::InputFloat("Azimuth##a.c_azimuth", &a.c_azimuth, 0, 0, "%.3f", disabledFlags);
		ImGui::InputFloat("Z##a.c_z", &a.c_z, 0, 0, "%.3f", disabledFlags);

		ImGui::End();
	}

	CrossProduct cp;
	void DrawCrossProduct() noexcept {
		ImGui::Begin("Cross Product");

		if (ImGui::InputFloat3("r ##cp.r", &cp.r[0]))
			cp.Calculate();

		DrawCopyPaste(&cp, &cp.r, COPY_BUTTON_STR(cp.r), PASTE_BUTTON_STR(cp.r));


		if (ImGui::InputFloat3("F##cp.F", &cp.F[0]))
			cp.Calculate();

		DrawCopyPaste(&cp, &cp.F, COPY_BUTTON_STR(cp.F), PASTE_BUTTON_STR(cp.F));

		ImGui::Text("");

		ImGui::InputFloat3("r x F##cp.M", &cp.M[0], "%.3f", disabledFlags);

		DrawCopy(&cp, &cp.M, COPY_BUTTON_STR(cp.M));

		ImGui::InputFloat("|r x F|##cp.Mmag", &cp.Mmag, 0, 0, "%.3f", disabledFlags);

		ImGui::End();
	}

	DotProduct dp;
	void DrawDotProduct() noexcept {
		ImGui::Begin("Dot Product");

		if (ImGui::InputFloat3("P##dp.P", &dp.P[0]))
			dp.Calculate();

		DrawCopyPaste(&dp, &dp.P, COPY_BUTTON_STR(dp.P), PASTE_BUTTON_STR(dp.P));

		if (ImGui::InputFloat3("Q##dp.Q", &dp.Q[0]))
			dp.Calculate();

		DrawCopyPaste(&dp, &dp.Q, COPY_BUTTON_STR(dp.Q), PASTE_BUTTON_STR(dp.Q));

		ImGui::Text("");

		ImGui::InputFloat("P dot Q##dp.R", &dp.R, 0, 0, "%.3f", disabledFlags);

		ImGui::End();
	}

	Distance d;
	void DrawDistance() noexcept {

		ImGui::Begin("Distance from Point to Line");

		if (ImGui::InputFloat3("P##d.P", &d.P[0]))
			d.Calculate();

		DrawCopyPaste(&d, &d.P, COPY_BUTTON_STR(d.P), PASTE_BUTTON_STR(d.P));

		ImGui::Text("");

		if (ImGui::InputFloat3("L0##d.L0", &d.L0[0]))
			d.Calculate();

		DrawCopyPaste(&d, &d.L0, COPY_BUTTON_STR(d.L0), PASTE_BUTTON_STR(d.L0));

		if (ImGui::InputFloat3("L   ##d.L", &d.L[0]))
			d.Calculate();

		DrawCopyPaste(&d, &d.L, COPY_BUTTON_STR(d.L), PASTE_BUTTON_STR(d.L));

		ImGui::Text("");

		ImGui::InputFloat("Distance##ddist", &d.dist, 0, 0, "%.3f", disabledFlags);

		ImGui::End();
	}

private:
	void DrawCopyPaste(Calculator* calculator, glm::vec3* var, const char* copystr, const char* pastestr) {
		ImGui::SameLine();
		if (ImGui::Button(copystr)) {
			clipboard = *var;
		}
		ImGui::SameLine();
		if (ImGui::Button(pastestr)) {
			*var = clipboard;
			calculator->Calculate();
		}
	}

	void DrawCopy(Calculator* calculator, glm::vec3* var, const char* copystr) {
		ImGui::SameLine();
		if (ImGui::Button(copystr)) {
			clipboard = *var;
		}
	}

	glm::vec3 clipboard;
	void DrawClipBoard() noexcept {
		ImGui::Begin("Clipboard");
		ImGui::InputFloat3("##Clipboard", &clipboard[0]);
		ImGui::End();
	}

public:
	virtual void OnUIRender() override
	{
		DrawClipBoard();
		DrawPositionVector();
		DrawUnitVector();
		DrawCrossProduct();
		DrawDotProduct();
		DrawDistance();
		DrawAngles();
	}
};


static int flag;

Walnut::Application* Walnut::CreateApplication(int argc, char** argv)
{
	Walnut::ApplicationSpecification spec;
	spec.Name = "Vector Calculator";

	Walnut::Application* app = new Walnut::Application(spec);
	app->PushLayer<VCLayer>();
	app->SetMenubarCallback([app]()
		{
			if (ImGui::BeginMenu("File"))
			{
				if (ImGui::MenuItem("Switch Theme")) {
					flag++;
					switch (flag % 3) {
					case 0:
						ImGui::StyleColorsDark();
						break;
					case 1:
						ImGui::StyleColorsLight();
						break;
					case 2:
						ImGui::StyleColorsClassic();
						break;
					}
				}

				if (ImGui::MenuItem("Exit"))
				{
					app->Close();
				}

				ImGui::EndMenu();
			}
		});
	return app;
}