/*
<%
setup_pybind11(cfg)
cfg['compiler_args'] = ['-std=c++11', '-I/usr/local/lib/python3.7/dist-packages/pybind11-2.4.3-py3.7.egg/']
cfg['dependencies'] = ['SimManager.hxx', 'DataManager.hxx', 'DataAnalyst.hxx', 'SimKernel.hxx', VortexCmplx.hxx', 'Complex.hxx']
cfg['sources'] = ['SimManager.cxx', 'DataManager.cxx', 'DataAnalyst.cxx', 'SimKernel.cxx', VortexCmplx.cxx', 'Complex.cxx']
%>
*/

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "SimManager.hxx"

namespace py = pybind11;

PYBIND11_MODULE(_vortex, m)
{
    py::class_<SimManager>(m, "SM")
    	.def(py::init<>())
        .def("setName", &SimManager::setName,
            py::arg("name"))
    	.def("addVtx", &SimManager::addVtx,
    		py::arg("x"),
    		py::arg("y"),
    		py::arg("circulation"),
    		py::arg("radius"),
            py::arg("fluid_id"))
    	.def("buildTimeSample", &SimManager::buildTimeSample,
    		py::arg("t0"),
    		py::arg("tEnd"),
    		py::arg("nbSteps"))
    	.def("setXPeriodicityTo", &SimManager::setXPeriodicityTo,
    		py::arg("x_periodicity"),
    		py::arg("x_period"))
    	.def("setYPeriodicityTo", &SimManager::setYPeriodicityTo,
    		py::arg("y_periodicity"),
    		py::arg("y_period"))
    	.def("chooseNumericalMethod", &SimManager::setMethodTo,
    		py::arg("method"))
    	.def("sim", &SimManager::sim,
            py::arg("nbThreads"))
        .def("getName", &SimManager::getName)
        .def("getNbVtx", &SimManager::getNbVtx)
        .def("getNbSteps", &SimManager::getNbSteps)
    	.def("getTimeList", &SimManager::getTimeVector)
    	.def("getXsAt", &SimManager::getXsAt,
    		py::arg("step"))
    	.def("getYsAt", &SimManager::getYsAt,
    		py::arg("step"))
    	.def("getUsAt", &SimManager::getUsAt,
    		py::arg("step"))
    	.def("getVsAt", &SimManager::getVsAt,
    		py::arg("step"))
        .def("getCirculationsAt", &SimManager::getCirculationsAt,
            py::arg("step"))
    	.def("getVtxXs", &SimManager::getVtxXs,
    		py::arg("vtxID"))
    	.def("getVtxYs", &SimManager::getVtxYs,
    		py::arg("vtxID"))
    	.def("getVtxUs", &SimManager::getVtxUs,
    		py::arg("vtxID"))
    	.def("getVtxVs", &SimManager::getVtxVs,
    		py::arg("vtxID"))
        .def("getVtxCirculations", &SimManager::getVtxCirculations,
            py::arg("vtxID"))
        .def("computeCompositionAt", &SimManager::computeCompositionAt,
            py::arg("x"),
            py::arg("y"),
            py::arg("step"),
            py::arg("radius"))
    	.def("computeVelocityAt", &SimManager::computeVelocityAt,
    		py::arg("x"),
    		py::arg("y"),
    		py::arg("step"),
    		py::arg("x_periodicity"),
    		py::arg("x_period"))
    	.def("computeVorticityAt", &SimManager::computeVorticityAt,
    		py::arg("x"),
    		py::arg("y"),
    		py::arg("step"),
    		py::arg("h"),
    		py::arg("x_periodicity"),
    		py::arg("x_period"))
        .def("computeHamiltonianEvolution", &SimManager::computeHamiltonianEvolution,
            py::arg("nbThreads"))
    	.def("saveSim", &SimManager::saveSim,
    		py::arg("fileNameWithPath"))
    	.def("loadSim", &SimManager::loadSim,
    		py::arg("fileNameWithPath"));
}