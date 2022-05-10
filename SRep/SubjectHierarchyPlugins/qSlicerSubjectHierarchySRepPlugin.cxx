/*==============================================================================

  Program: 3D Slicer

  Copyright (c) Laboratory for Percutaneous Surgery (PerkLab)
  Queen's University, Kingston, ON, Canada. All Rights Reserved.

  See COPYRIGHT.txt
  or http://www.slicer.org/copyright/copyright.txt for details.

  Unless required by applicable law or agreed to in writing, software
  distributed under the License is distributed on an "AS IS" BASIS,
  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
  See the License for the specific language governing permissions and
  limitations under the License.

  This file was originally developed by Csaba Pinter, PerkLab, Queen's University
  and was supported through the Applied Cancer Research Unit program of Cancer Care
  Ontario with funds provided by the Ontario Ministry of Health and Long-Term Care

==============================================================================*/

// SubjectHierarchy Plugins includes
#include "qSlicerSubjectHierarchyDefaultPlugin.h"
#include "qSlicerSubjectHierarchyFolderPlugin.h"
#include "qSlicerSubjectHierarchySRepPlugin.h"
#include "qSlicerSubjectHierarchyPluginHandler.h"

// Terminologies includes
#include "qSlicerTerminologyItemDelegate.h"
#include "vtkSlicerTerminologiesModuleLogic.h"

// MRML includes
#include <vtkMRMLScene.h>
#include <vtkMRMLSRepNode.h>
#include <vtkMRMLSRepDisplayNode.h>

// vtkSegmentationCore includes
#include <vtkSegment.h>

// VTK includes
#include <vtkObjectFactory.h>
#include <vtkSmartPointer.h>
#include <vtkPolyData.h>

// Qt includes
#include <QDebug>
#include <QStandardItem>
#include <QAction>

// Slicer includes
#include "qSlicerAbstractModuleWidget.h"

//-----------------------------------------------------------------------------
/// \ingroup Slicer_QtModules_SubjectHierarchy_Plugins
class qSlicerSubjectHierarchySRepPluginPrivate: public QObject
{
  Q_DECLARE_PUBLIC(qSlicerSubjectHierarchySRepPlugin);
protected:
  qSlicerSubjectHierarchySRepPlugin* const q_ptr;
public:
  qSlicerSubjectHierarchySRepPluginPrivate(qSlicerSubjectHierarchySRepPlugin& object);
  ~qSlicerSubjectHierarchySRepPluginPrivate() override;
  void init();
public:
  QIcon SRepIcon;
};

//-----------------------------------------------------------------------------
// qSlicerSubjectHierarchySRepPluginPrivate methods

//-----------------------------------------------------------------------------
qSlicerSubjectHierarchySRepPluginPrivate::qSlicerSubjectHierarchySRepPluginPrivate(qSlicerSubjectHierarchySRepPlugin& object)
: q_ptr(&object)
, SRepIcon(QIcon(":Icons/SRep.png"))
{
}

//------------------------------------------------------------------------------
void qSlicerSubjectHierarchySRepPluginPrivate::init()
{
}

//-----------------------------------------------------------------------------
qSlicerSubjectHierarchySRepPluginPrivate::~qSlicerSubjectHierarchySRepPluginPrivate() = default;

//-----------------------------------------------------------------------------
// qSlicerSubjectHierarchySRepPlugin methods

//-----------------------------------------------------------------------------
qSlicerSubjectHierarchySRepPlugin::qSlicerSubjectHierarchySRepPlugin(QObject* parent)
 : Superclass(parent)
 , d_ptr( new qSlicerSubjectHierarchySRepPluginPrivate(*this) )
{
  this->m_Name = QString("SRep");

  Q_D(qSlicerSubjectHierarchySRepPlugin);
  d->init();
}

//-----------------------------------------------------------------------------
qSlicerSubjectHierarchySRepPlugin::~qSlicerSubjectHierarchySRepPlugin() = default;

//----------------------------------------------------------------------------
double qSlicerSubjectHierarchySRepPlugin::canAddNodeToSubjectHierarchy(
  vtkMRMLNode* node, vtkIdType parentItemID/*=vtkMRMLSubjectHierarchyNode::INVALID_ITEM_ID*/)const
{
  Q_UNUSED(parentItemID);
  if (!node) {
    qCritical() << Q_FUNC_INFO << ": Input node is NULL";
    return 0.0;
  }
  else if (node->IsA("vtkMRMLSRepNode")) {
    // Node is a srep
    return 0.5;
  }
  return 0.0;
}

//---------------------------------------------------------------------------
double qSlicerSubjectHierarchySRepPlugin::canOwnSubjectHierarchyItem(vtkIdType itemID)const
{
  if (itemID == vtkMRMLSubjectHierarchyNode::INVALID_ITEM_ID) {
    qCritical() << Q_FUNC_INFO << ": Invalid input item";
    return 0.0;
  }
  vtkMRMLSubjectHierarchyNode* shNode = qSlicerSubjectHierarchyPluginHandler::instance()->subjectHierarchyNode();
  if (!shNode) {
    qCritical() << Q_FUNC_INFO << ": Failed to access subject hierarchy node";
    return 0.0;
  }

  // SRep
  vtkMRMLNode* associatedNode = shNode->GetItemDataNode(itemID);
  if (associatedNode && associatedNode->IsA("vtkMRMLSRepNode")) {
    return 0.5; // There may be other plugins that can handle special srep better
  }

  return 0.0;
}

//---------------------------------------------------------------------------
const QString qSlicerSubjectHierarchySRepPlugin::roleForPlugin()const
{
  return "SRep";
}

//---------------------------------------------------------------------------
QIcon qSlicerSubjectHierarchySRepPlugin::icon(vtkIdType itemID)
{
  Q_D(qSlicerSubjectHierarchySRepPlugin);

  if (itemID == vtkMRMLSubjectHierarchyNode::INVALID_ITEM_ID) {
    qCritical() << Q_FUNC_INFO << ": Invalid input item";
    return QIcon();
  }

  if (this->canOwnSubjectHierarchyItem(itemID)) {
    return d->SRepIcon;
  }

  // Item unknown by plugin
  return QIcon();
}

//---------------------------------------------------------------------------
QIcon qSlicerSubjectHierarchySRepPlugin::visibilityIcon(int visible)
{
  // Have the default plugin (which is not registered) take care of this
  return qSlicerSubjectHierarchyPluginHandler::instance()->defaultPlugin()->visibilityIcon(visible);
}

//-----------------------------------------------------------------------------
QString qSlicerSubjectHierarchySRepPlugin::tooltip(vtkIdType itemID)const
{
  if (itemID == vtkMRMLSubjectHierarchyNode::INVALID_ITEM_ID) {
    qCritical() << Q_FUNC_INFO << ": Invalid input item";
    return QString("Invalid");
  }
  vtkMRMLSubjectHierarchyNode* shNode = qSlicerSubjectHierarchyPluginHandler::instance()->subjectHierarchyNode();
  if (!shNode) {
    qCritical() << Q_FUNC_INFO << ": Failed to access subject hierarchy node";
    return QString("Error");
  }

  // Get basic tooltip from abstract plugin
  QString tooltipString = Superclass::tooltip(itemID);

  vtkMRMLSRepNode* srepNode = vtkMRMLSRepNode::SafeDownCast(shNode->GetItemDataNode(itemID));

  auto displayNode = srepNode ? vtkMRMLSRepDisplayNode::SafeDownCast(srepNode->GetDisplayNode()) : nullptr;
  if (srepNode && displayNode) {
    tooltipString.append( QString(" (Has SRep: %1)")
      .arg((srepNode->GetSRep() && !srepNode->GetSRep()->IsEmpty()) ? "YES" : "NO") );
  }
  else {
    tooltipString.append(" !Invalid srep");
  }

  return tooltipString;
}

//-----------------------------------------------------------------------------
void qSlicerSubjectHierarchySRepPlugin::setDisplayColor(vtkIdType itemID, QColor color, QMap<int, QVariant> terminologyMetaData)
{
  if (itemID == vtkMRMLSubjectHierarchyNode::INVALID_ITEM_ID) {
    qCritical() << Q_FUNC_INFO << ": Invalid input item";
    return;
  }
  vtkMRMLSubjectHierarchyNode* shNode = qSlicerSubjectHierarchyPluginHandler::instance()->subjectHierarchyNode();
  if (!shNode) {
    qCritical() << Q_FUNC_INFO << ": Failed to access subject hierarchy node";
    return;
  }

  // Get srep node and display node
  vtkMRMLSRepNode* srepNode = vtkMRMLSRepNode::SafeDownCast(shNode->GetItemDataNode(itemID));
  if (!srepNode) {
    qCritical() << Q_FUNC_INFO << ": Unable to find srep node for subject hierarchy item " << shNode->GetItemName(itemID).c_str();
    return;
  }
  vtkMRMLSRepDisplayNode* displayNode = vtkMRMLSRepDisplayNode::SafeDownCast(srepNode->GetDisplayNode());
  if (!displayNode) {
    qCritical() << Q_FUNC_INFO << ": No display node for srep";
    return;
  }

  // Set terminology metadata
  if (terminologyMetaData.contains(qSlicerTerminologyItemDelegate::TerminologyRole)) {
    srepNode->SetAttribute(vtkSegment::GetTerminologyEntryTagName(),
      terminologyMetaData[qSlicerTerminologyItemDelegate::TerminologyRole].toString().toUtf8().constData() );
  }
  if (terminologyMetaData.contains(qSlicerTerminologyItemDelegate::NameRole)) {
    srepNode->SetName(
      terminologyMetaData[qSlicerTerminologyItemDelegate::NameRole].toString().toUtf8().constData() );
  }
  if (terminologyMetaData.contains(qSlicerTerminologyItemDelegate::NameAutoGeneratedRole)) {
    srepNode->SetAttribute( vtkSlicerTerminologiesModuleLogic::GetNameAutoGeneratedAttributeName(),
      terminologyMetaData[qSlicerTerminologyItemDelegate::NameAutoGeneratedRole].toString().toUtf8().constData() );
  }
  if (terminologyMetaData.contains(qSlicerTerminologyItemDelegate::ColorAutoGeneratedRole)) {
    srepNode->SetAttribute( vtkSlicerTerminologiesModuleLogic::GetColorAutoGeneratedAttributeName(),
      terminologyMetaData[qSlicerTerminologyItemDelegate::ColorAutoGeneratedRole].toString().toUtf8().constData() );
  }

  // Set color
  double* oldColorArray = displayNode->GetColor();
  QColor oldColor = QColor::fromRgbF(oldColorArray[0], oldColorArray[1], oldColorArray[2]);
  if (oldColor != color) {
    displayNode->SetColor(color.redF(), color.greenF(), color.blueF());

    // Trigger update of color swatch
    shNode->ItemModified(itemID);
  }
}

//-----------------------------------------------------------------------------
QColor qSlicerSubjectHierarchySRepPlugin::getDisplayColor(vtkIdType itemID, QMap<int, QVariant> &terminologyMetaData)const
{
  if (itemID == vtkMRMLSubjectHierarchyNode::INVALID_ITEM_ID) {
    qCritical() << Q_FUNC_INFO << ": Invalid input item";
    return QColor(0,0,0,0);
  }
  vtkMRMLSubjectHierarchyNode* shNode = qSlicerSubjectHierarchyPluginHandler::instance()->subjectHierarchyNode();
  if (!shNode) {
    qCritical() << Q_FUNC_INFO << ": Failed to access subject hierarchy node";
    return QColor(0,0,0,0);
  }
  vtkMRMLScene* scene = qSlicerSubjectHierarchyPluginHandler::instance()->mrmlScene();
  if (!scene) {
    qCritical() << Q_FUNC_INFO << ": Invalid MRML scene";
    return QColor(0,0,0,0);
  }

  if (scene->IsImporting()) {
    // During import SH node may be created before the segmentation is read into the scene,
    // so don't attempt to access the segment yet
    return QColor(0,0,0,0);
  }

  // Get srep node and display node
  vtkMRMLSRepNode* srepNode = vtkMRMLSRepNode::SafeDownCast(shNode->GetItemDataNode(itemID));
  if (!srepNode) {
    qCritical() << Q_FUNC_INFO << ": Unable to find srep node for subject hierarchy item " << shNode->GetItemName(itemID).c_str();
    return QColor(0,0,0,0);
  }
  vtkMRMLSRepDisplayNode* displayNode = vtkMRMLSRepDisplayNode::SafeDownCast(srepNode->GetDisplayNode());
  if (!displayNode) {
    return QColor(0,0,0,0);
  }

  // Get terminology metadata
  terminologyMetaData.clear();
  terminologyMetaData[qSlicerTerminologyItemDelegate::TerminologyRole] =
    srepNode->GetAttribute(vtkSegment::GetTerminologyEntryTagName());
  terminologyMetaData[qSlicerTerminologyItemDelegate::NameRole] = srepNode->GetName();
  // If auto generated flags are not initialized, then set them to the default
  // (color: on, name: off - this way color will be set from the selector but name will not)
  bool nameAutoGenerated = false;
  if (srepNode->GetAttribute(vtkSlicerTerminologiesModuleLogic::GetNameAutoGeneratedAttributeName())) {
    nameAutoGenerated = QVariant(srepNode->GetAttribute(vtkSlicerTerminologiesModuleLogic::GetNameAutoGeneratedAttributeName())).toBool();
  }
  terminologyMetaData[qSlicerTerminologyItemDelegate::NameAutoGeneratedRole] = nameAutoGenerated;
  bool colorAutoGenerated = true;
  if (srepNode->GetAttribute(vtkSlicerTerminologiesModuleLogic::GetColorAutoGeneratedAttributeName())) {
    colorAutoGenerated = QVariant(srepNode->GetAttribute(vtkSlicerTerminologiesModuleLogic::GetColorAutoGeneratedAttributeName())).toBool();
  }
  terminologyMetaData[qSlicerTerminologyItemDelegate::ColorAutoGeneratedRole] = colorAutoGenerated;

  // Get and return color
  double* colorArray = displayNode->GetColor();
  return QColor::fromRgbF(colorArray[0], colorArray[1], colorArray[2]);
}
