#-----------------------------------------------------------------------------
set(MODULE_NAME SRepVisualizer)

#-----------------------------------------------------------------------------
set(MODULE_PYTHON_SCRIPTS
  ${MODULE_NAME}.py
  LegacyTransformer/__init__.py
  LegacyTransformer/srep.py
  LegacyTransformer/srep_io.py
  LegacyTransformer/legacyTransformer.py
  LegacyTransformer/ElementTree_pretty.py
  )

set(MODULE_PYTHON_RESOURCES
  Resources/Icons/${MODULE_NAME}.png
  )

#-----------------------------------------------------------------------------
slicerMacroBuildScriptedModule(
  NAME ${MODULE_NAME}
  SCRIPTS ${MODULE_PYTHON_SCRIPTS}
  RESOURCES ${MODULE_PYTHON_RESOURCES}
  WITH_GENERIC_TESTS
  )

#-----------------------------------------------------------------------------
if(BUILD_TESTING)

  # Additional build-time testing
  add_subdirectory(Testing)
endif()
