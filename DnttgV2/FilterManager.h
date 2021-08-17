#pragma once
#include "pandaFramework.h"
#include "pandaSystem.h"

class FilterManager
{
public:
	FilterManager(WindowFramework *window, NodePath c);
private:
	WindowFramework * win;
	NodePath cam;
};

