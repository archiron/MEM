// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME tmpdIslc6_amd64_gcc530dIsrcdIMEMDataFormatsdIMEMDataFormatsdIsrcdIMEMDataFormatsMEMDataFormatsdIadIMEMDataFormatsMEMDataFormats_xr

/*******************************************************************/
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#define G__DICTIONARY
#include "RConfig.h"
#include "TClass.h"
#include "TDictAttributeMap.h"
#include "TInterpreter.h"
#include "TROOT.h"
#include "TBuffer.h"
#include "TMemberInspector.h"
#include "TInterpreter.h"
#include "TVirtualMutex.h"
#include "TError.h"

#ifndef G__ROOT
#define G__ROOT
#endif

#include "RtypesImp.h"
#include "TIsAProxy.h"
#include "TFileMergeInfo.h"
#include <algorithm>
#include "TCollectionProxyInfo.h"
/*******************************************************************/

#include "TDataMember.h"

// Since CINT ignores the std namespace, we need to do so in this file.
namespace std {} using namespace std;

// Header files passed as explicit arguments
#include "src/MEMDataFormats/MEMDataFormats/src/classes.h"

// Header files passed via #pragma extra_include

namespace ROOT {
   static TClass *recocLcLMEMResult_Dictionary();
   static void recocLcLMEMResult_TClassManip(TClass*);
   static void *new_recocLcLMEMResult(void *p = 0);
   static void *newArray_recocLcLMEMResult(Long_t size, void *p);
   static void delete_recocLcLMEMResult(void *p);
   static void deleteArray_recocLcLMEMResult(void *p);
   static void destruct_recocLcLMEMResult(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::reco::MEMResult*)
   {
      ::reco::MEMResult *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::reco::MEMResult));
      static ::ROOT::TGenericClassInfo 
         instance("reco::MEMResult", 15, "MEMDataFormats/MEMDataFormats/interface/MEMResult.h", 69,
                  typeid(::reco::MEMResult), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &recocLcLMEMResult_Dictionary, isa_proxy, 8,
                  sizeof(::reco::MEMResult) );
      instance.SetNew(&new_recocLcLMEMResult);
      instance.SetNewArray(&newArray_recocLcLMEMResult);
      instance.SetDelete(&delete_recocLcLMEMResult);
      instance.SetDeleteArray(&deleteArray_recocLcLMEMResult);
      instance.SetDestructor(&destruct_recocLcLMEMResult);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::reco::MEMResult*)
   {
      return GenerateInitInstanceLocal((::reco::MEMResult*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::reco::MEMResult*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *recocLcLMEMResult_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::reco::MEMResult*)0x0)->GetClass();
      recocLcLMEMResult_TClassManip(theClass);
   return theClass;
   }

   static void recocLcLMEMResult_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *edmcLcLWrapperlEvectorlErecocLcLMEMResultgRsPgR_Dictionary();
   static void edmcLcLWrapperlEvectorlErecocLcLMEMResultgRsPgR_TClassManip(TClass*);
   static void *new_edmcLcLWrapperlEvectorlErecocLcLMEMResultgRsPgR(void *p = 0);
   static void *newArray_edmcLcLWrapperlEvectorlErecocLcLMEMResultgRsPgR(Long_t size, void *p);
   static void delete_edmcLcLWrapperlEvectorlErecocLcLMEMResultgRsPgR(void *p);
   static void deleteArray_edmcLcLWrapperlEvectorlErecocLcLMEMResultgRsPgR(void *p);
   static void destruct_edmcLcLWrapperlEvectorlErecocLcLMEMResultgRsPgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::edm::Wrapper<vector<reco::MEMResult> >*)
   {
      ::edm::Wrapper<vector<reco::MEMResult> > *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::edm::Wrapper<vector<reco::MEMResult> >));
      static ::ROOT::TGenericClassInfo 
         instance("edm::Wrapper<vector<reco::MEMResult> >", ::edm::Wrapper<vector<reco::MEMResult> >::Class_Version(), "DataFormats/Common/interface/Wrapper.h", 25,
                  typeid(::edm::Wrapper<vector<reco::MEMResult> >), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &edmcLcLWrapperlEvectorlErecocLcLMEMResultgRsPgR_Dictionary, isa_proxy, 0,
                  sizeof(::edm::Wrapper<vector<reco::MEMResult> >) );
      instance.SetNew(&new_edmcLcLWrapperlEvectorlErecocLcLMEMResultgRsPgR);
      instance.SetNewArray(&newArray_edmcLcLWrapperlEvectorlErecocLcLMEMResultgRsPgR);
      instance.SetDelete(&delete_edmcLcLWrapperlEvectorlErecocLcLMEMResultgRsPgR);
      instance.SetDeleteArray(&deleteArray_edmcLcLWrapperlEvectorlErecocLcLMEMResultgRsPgR);
      instance.SetDestructor(&destruct_edmcLcLWrapperlEvectorlErecocLcLMEMResultgRsPgR);

      ::ROOT::AddClassAlternate("edm::Wrapper<vector<reco::MEMResult> >","edm::Wrapper<std::vector<reco::MEMResult> >");
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::edm::Wrapper<vector<reco::MEMResult> >*)
   {
      return GenerateInitInstanceLocal((::edm::Wrapper<vector<reco::MEMResult> >*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::edm::Wrapper<vector<reco::MEMResult> >*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *edmcLcLWrapperlEvectorlErecocLcLMEMResultgRsPgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::edm::Wrapper<vector<reco::MEMResult> >*)0x0)->GetClass();
      edmcLcLWrapperlEvectorlErecocLcLMEMResultgRsPgR_TClassManip(theClass);
   return theClass;
   }

   static void edmcLcLWrapperlEvectorlErecocLcLMEMResultgRsPgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *edmcLcLReflEvectorlErecocLcLMEMResultgRcOrecocLcLMEMResultcOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMEMResultgRcOrecocLcLMEMResultgRsPgR_Dictionary();
   static void edmcLcLReflEvectorlErecocLcLMEMResultgRcOrecocLcLMEMResultcOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMEMResultgRcOrecocLcLMEMResultgRsPgR_TClassManip(TClass*);
   static void *new_edmcLcLReflEvectorlErecocLcLMEMResultgRcOrecocLcLMEMResultcOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMEMResultgRcOrecocLcLMEMResultgRsPgR(void *p = 0);
   static void *newArray_edmcLcLReflEvectorlErecocLcLMEMResultgRcOrecocLcLMEMResultcOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMEMResultgRcOrecocLcLMEMResultgRsPgR(Long_t size, void *p);
   static void delete_edmcLcLReflEvectorlErecocLcLMEMResultgRcOrecocLcLMEMResultcOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMEMResultgRcOrecocLcLMEMResultgRsPgR(void *p);
   static void deleteArray_edmcLcLReflEvectorlErecocLcLMEMResultgRcOrecocLcLMEMResultcOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMEMResultgRcOrecocLcLMEMResultgRsPgR(void *p);
   static void destruct_edmcLcLReflEvectorlErecocLcLMEMResultgRcOrecocLcLMEMResultcOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMEMResultgRcOrecocLcLMEMResultgRsPgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::edm::Ref<vector<reco::MEMResult>,reco::MEMResult,edm::refhelper::FindUsingAdvance<vector<reco::MEMResult>,reco::MEMResult> >*)
   {
      ::edm::Ref<vector<reco::MEMResult>,reco::MEMResult,edm::refhelper::FindUsingAdvance<vector<reco::MEMResult>,reco::MEMResult> > *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::edm::Ref<vector<reco::MEMResult>,reco::MEMResult,edm::refhelper::FindUsingAdvance<vector<reco::MEMResult>,reco::MEMResult> >));
      static ::ROOT::TGenericClassInfo 
         instance("edm::Ref<vector<reco::MEMResult>,reco::MEMResult,edm::refhelper::FindUsingAdvance<vector<reco::MEMResult>,reco::MEMResult> >", ::edm::Ref<vector<reco::MEMResult>,reco::MEMResult,edm::refhelper::FindUsingAdvance<vector<reco::MEMResult>,reco::MEMResult> >::Class_Version(), "DataFormats/Common/interface/Ref.h", 305,
                  typeid(::edm::Ref<vector<reco::MEMResult>,reco::MEMResult,edm::refhelper::FindUsingAdvance<vector<reco::MEMResult>,reco::MEMResult> >), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &edmcLcLReflEvectorlErecocLcLMEMResultgRcOrecocLcLMEMResultcOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMEMResultgRcOrecocLcLMEMResultgRsPgR_Dictionary, isa_proxy, 0,
                  sizeof(::edm::Ref<vector<reco::MEMResult>,reco::MEMResult,edm::refhelper::FindUsingAdvance<vector<reco::MEMResult>,reco::MEMResult> >) );
      instance.SetNew(&new_edmcLcLReflEvectorlErecocLcLMEMResultgRcOrecocLcLMEMResultcOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMEMResultgRcOrecocLcLMEMResultgRsPgR);
      instance.SetNewArray(&newArray_edmcLcLReflEvectorlErecocLcLMEMResultgRcOrecocLcLMEMResultcOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMEMResultgRcOrecocLcLMEMResultgRsPgR);
      instance.SetDelete(&delete_edmcLcLReflEvectorlErecocLcLMEMResultgRcOrecocLcLMEMResultcOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMEMResultgRcOrecocLcLMEMResultgRsPgR);
      instance.SetDeleteArray(&deleteArray_edmcLcLReflEvectorlErecocLcLMEMResultgRcOrecocLcLMEMResultcOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMEMResultgRcOrecocLcLMEMResultgRsPgR);
      instance.SetDestructor(&destruct_edmcLcLReflEvectorlErecocLcLMEMResultgRcOrecocLcLMEMResultcOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMEMResultgRcOrecocLcLMEMResultgRsPgR);

      ::ROOT::AddClassAlternate("edm::Ref<vector<reco::MEMResult>,reco::MEMResult,edm::refhelper::FindUsingAdvance<vector<reco::MEMResult>,reco::MEMResult> >","edm::Ref<std::vector<reco::MEMResult>,reco::MEMResult,edm::refhelper::FindUsingAdvance<std::vector<reco::MEMResult>,reco::MEMResult> >");
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::edm::Ref<vector<reco::MEMResult>,reco::MEMResult,edm::refhelper::FindUsingAdvance<vector<reco::MEMResult>,reco::MEMResult> >*)
   {
      return GenerateInitInstanceLocal((::edm::Ref<vector<reco::MEMResult>,reco::MEMResult,edm::refhelper::FindUsingAdvance<vector<reco::MEMResult>,reco::MEMResult> >*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::edm::Ref<vector<reco::MEMResult>,reco::MEMResult,edm::refhelper::FindUsingAdvance<vector<reco::MEMResult>,reco::MEMResult> >*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *edmcLcLReflEvectorlErecocLcLMEMResultgRcOrecocLcLMEMResultcOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMEMResultgRcOrecocLcLMEMResultgRsPgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::edm::Ref<vector<reco::MEMResult>,reco::MEMResult,edm::refhelper::FindUsingAdvance<vector<reco::MEMResult>,reco::MEMResult> >*)0x0)->GetClass();
      edmcLcLReflEvectorlErecocLcLMEMResultgRcOrecocLcLMEMResultcOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMEMResultgRcOrecocLcLMEMResultgRsPgR_TClassManip(theClass);
   return theClass;
   }

   static void edmcLcLReflEvectorlErecocLcLMEMResultgRcOrecocLcLMEMResultcOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMEMResultgRcOrecocLcLMEMResultgRsPgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *edmcLcLRefProdlEvectorlErecocLcLMEMResultgRsPgR_Dictionary();
   static void edmcLcLRefProdlEvectorlErecocLcLMEMResultgRsPgR_TClassManip(TClass*);
   static void *new_edmcLcLRefProdlEvectorlErecocLcLMEMResultgRsPgR(void *p = 0);
   static void *newArray_edmcLcLRefProdlEvectorlErecocLcLMEMResultgRsPgR(Long_t size, void *p);
   static void delete_edmcLcLRefProdlEvectorlErecocLcLMEMResultgRsPgR(void *p);
   static void deleteArray_edmcLcLRefProdlEvectorlErecocLcLMEMResultgRsPgR(void *p);
   static void destruct_edmcLcLRefProdlEvectorlErecocLcLMEMResultgRsPgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::edm::RefProd<vector<reco::MEMResult> >*)
   {
      ::edm::RefProd<vector<reco::MEMResult> > *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::edm::RefProd<vector<reco::MEMResult> >));
      static ::ROOT::TGenericClassInfo 
         instance("edm::RefProd<vector<reco::MEMResult> >", ::edm::RefProd<vector<reco::MEMResult> >::Class_Version(), "DataFormats/Common/interface/RefProd.h", 55,
                  typeid(::edm::RefProd<vector<reco::MEMResult> >), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &edmcLcLRefProdlEvectorlErecocLcLMEMResultgRsPgR_Dictionary, isa_proxy, 0,
                  sizeof(::edm::RefProd<vector<reco::MEMResult> >) );
      instance.SetNew(&new_edmcLcLRefProdlEvectorlErecocLcLMEMResultgRsPgR);
      instance.SetNewArray(&newArray_edmcLcLRefProdlEvectorlErecocLcLMEMResultgRsPgR);
      instance.SetDelete(&delete_edmcLcLRefProdlEvectorlErecocLcLMEMResultgRsPgR);
      instance.SetDeleteArray(&deleteArray_edmcLcLRefProdlEvectorlErecocLcLMEMResultgRsPgR);
      instance.SetDestructor(&destruct_edmcLcLRefProdlEvectorlErecocLcLMEMResultgRsPgR);

      ::ROOT::AddClassAlternate("edm::RefProd<vector<reco::MEMResult> >","edm::RefProd<std::vector<reco::MEMResult> >");
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::edm::RefProd<vector<reco::MEMResult> >*)
   {
      return GenerateInitInstanceLocal((::edm::RefProd<vector<reco::MEMResult> >*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::edm::RefProd<vector<reco::MEMResult> >*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *edmcLcLRefProdlEvectorlErecocLcLMEMResultgRsPgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::edm::RefProd<vector<reco::MEMResult> >*)0x0)->GetClass();
      edmcLcLRefProdlEvectorlErecocLcLMEMResultgRsPgR_TClassManip(theClass);
   return theClass;
   }

   static void edmcLcLRefProdlEvectorlErecocLcLMEMResultgRsPgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *edmcLcLRefVectorlEvectorlErecocLcLMEMResultgRcOrecocLcLMEMResultcOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMEMResultgRcOrecocLcLMEMResultgRsPgR_Dictionary();
   static void edmcLcLRefVectorlEvectorlErecocLcLMEMResultgRcOrecocLcLMEMResultcOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMEMResultgRcOrecocLcLMEMResultgRsPgR_TClassManip(TClass*);
   static void *new_edmcLcLRefVectorlEvectorlErecocLcLMEMResultgRcOrecocLcLMEMResultcOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMEMResultgRcOrecocLcLMEMResultgRsPgR(void *p = 0);
   static void *newArray_edmcLcLRefVectorlEvectorlErecocLcLMEMResultgRcOrecocLcLMEMResultcOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMEMResultgRcOrecocLcLMEMResultgRsPgR(Long_t size, void *p);
   static void delete_edmcLcLRefVectorlEvectorlErecocLcLMEMResultgRcOrecocLcLMEMResultcOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMEMResultgRcOrecocLcLMEMResultgRsPgR(void *p);
   static void deleteArray_edmcLcLRefVectorlEvectorlErecocLcLMEMResultgRcOrecocLcLMEMResultcOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMEMResultgRcOrecocLcLMEMResultgRsPgR(void *p);
   static void destruct_edmcLcLRefVectorlEvectorlErecocLcLMEMResultgRcOrecocLcLMEMResultcOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMEMResultgRcOrecocLcLMEMResultgRsPgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::edm::RefVector<vector<reco::MEMResult>,reco::MEMResult,edm::refhelper::FindUsingAdvance<vector<reco::MEMResult>,reco::MEMResult> >*)
   {
      ::edm::RefVector<vector<reco::MEMResult>,reco::MEMResult,edm::refhelper::FindUsingAdvance<vector<reco::MEMResult>,reco::MEMResult> > *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::edm::RefVector<vector<reco::MEMResult>,reco::MEMResult,edm::refhelper::FindUsingAdvance<vector<reco::MEMResult>,reco::MEMResult> >));
      static ::ROOT::TGenericClassInfo 
         instance("edm::RefVector<vector<reco::MEMResult>,reco::MEMResult,edm::refhelper::FindUsingAdvance<vector<reco::MEMResult>,reco::MEMResult> >", ::edm::RefVector<vector<reco::MEMResult>,reco::MEMResult,edm::refhelper::FindUsingAdvance<vector<reco::MEMResult>,reco::MEMResult> >::Class_Version(), "DataFormats/Common/interface/RefVector.h", 33,
                  typeid(::edm::RefVector<vector<reco::MEMResult>,reco::MEMResult,edm::refhelper::FindUsingAdvance<vector<reco::MEMResult>,reco::MEMResult> >), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &edmcLcLRefVectorlEvectorlErecocLcLMEMResultgRcOrecocLcLMEMResultcOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMEMResultgRcOrecocLcLMEMResultgRsPgR_Dictionary, isa_proxy, 0,
                  sizeof(::edm::RefVector<vector<reco::MEMResult>,reco::MEMResult,edm::refhelper::FindUsingAdvance<vector<reco::MEMResult>,reco::MEMResult> >) );
      instance.SetNew(&new_edmcLcLRefVectorlEvectorlErecocLcLMEMResultgRcOrecocLcLMEMResultcOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMEMResultgRcOrecocLcLMEMResultgRsPgR);
      instance.SetNewArray(&newArray_edmcLcLRefVectorlEvectorlErecocLcLMEMResultgRcOrecocLcLMEMResultcOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMEMResultgRcOrecocLcLMEMResultgRsPgR);
      instance.SetDelete(&delete_edmcLcLRefVectorlEvectorlErecocLcLMEMResultgRcOrecocLcLMEMResultcOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMEMResultgRcOrecocLcLMEMResultgRsPgR);
      instance.SetDeleteArray(&deleteArray_edmcLcLRefVectorlEvectorlErecocLcLMEMResultgRcOrecocLcLMEMResultcOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMEMResultgRcOrecocLcLMEMResultgRsPgR);
      instance.SetDestructor(&destruct_edmcLcLRefVectorlEvectorlErecocLcLMEMResultgRcOrecocLcLMEMResultcOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMEMResultgRcOrecocLcLMEMResultgRsPgR);

      ::ROOT::AddClassAlternate("edm::RefVector<vector<reco::MEMResult>,reco::MEMResult,edm::refhelper::FindUsingAdvance<vector<reco::MEMResult>,reco::MEMResult> >","edm::RefVector<std::vector<reco::MEMResult>,reco::MEMResult,edm::refhelper::FindUsingAdvance<std::vector<reco::MEMResult>,reco::MEMResult> >");
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::edm::RefVector<vector<reco::MEMResult>,reco::MEMResult,edm::refhelper::FindUsingAdvance<vector<reco::MEMResult>,reco::MEMResult> >*)
   {
      return GenerateInitInstanceLocal((::edm::RefVector<vector<reco::MEMResult>,reco::MEMResult,edm::refhelper::FindUsingAdvance<vector<reco::MEMResult>,reco::MEMResult> >*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::edm::RefVector<vector<reco::MEMResult>,reco::MEMResult,edm::refhelper::FindUsingAdvance<vector<reco::MEMResult>,reco::MEMResult> >*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *edmcLcLRefVectorlEvectorlErecocLcLMEMResultgRcOrecocLcLMEMResultcOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMEMResultgRcOrecocLcLMEMResultgRsPgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::edm::RefVector<vector<reco::MEMResult>,reco::MEMResult,edm::refhelper::FindUsingAdvance<vector<reco::MEMResult>,reco::MEMResult> >*)0x0)->GetClass();
      edmcLcLRefVectorlEvectorlErecocLcLMEMResultgRcOrecocLcLMEMResultcOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMEMResultgRcOrecocLcLMEMResultgRsPgR_TClassManip(theClass);
   return theClass;
   }

   static void edmcLcLRefVectorlEvectorlErecocLcLMEMResultgRcOrecocLcLMEMResultcOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMEMResultgRcOrecocLcLMEMResultgRsPgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *recocLcLIntegralResult_Dictionary();
   static void recocLcLIntegralResult_TClassManip(TClass*);
   static void *new_recocLcLIntegralResult(void *p = 0);
   static void *newArray_recocLcLIntegralResult(Long_t size, void *p);
   static void delete_recocLcLIntegralResult(void *p);
   static void deleteArray_recocLcLIntegralResult(void *p);
   static void destruct_recocLcLIntegralResult(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::reco::IntegralResult*)
   {
      ::reco::IntegralResult *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::reco::IntegralResult));
      static ::ROOT::TGenericClassInfo 
         instance("reco::IntegralResult", "MEMDataFormats/MEMDataFormats/interface/IntegralResult.h", 62,
                  typeid(::reco::IntegralResult), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &recocLcLIntegralResult_Dictionary, isa_proxy, 0,
                  sizeof(::reco::IntegralResult) );
      instance.SetNew(&new_recocLcLIntegralResult);
      instance.SetNewArray(&newArray_recocLcLIntegralResult);
      instance.SetDelete(&delete_recocLcLIntegralResult);
      instance.SetDeleteArray(&deleteArray_recocLcLIntegralResult);
      instance.SetDestructor(&destruct_recocLcLIntegralResult);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::reco::IntegralResult*)
   {
      return GenerateInitInstanceLocal((::reco::IntegralResult*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::reco::IntegralResult*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *recocLcLIntegralResult_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::reco::IntegralResult*)0x0)->GetClass();
      recocLcLIntegralResult_TClassManip(theClass);
   return theClass;
   }

   static void recocLcLIntegralResult_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *RunConfig_Dictionary();
   static void RunConfig_TClassManip(TClass*);
   static void *new_RunConfig(void *p = 0);
   static void *newArray_RunConfig(Long_t size, void *p);
   static void delete_RunConfig(void *p);
   static void deleteArray_RunConfig(void *p);
   static void destruct_RunConfig(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::RunConfig*)
   {
      ::RunConfig *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::RunConfig));
      static ::ROOT::TGenericClassInfo 
         instance("RunConfig", "MEMDataFormats/MEMDataFormats/interface/RunConfig.h", 28,
                  typeid(::RunConfig), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &RunConfig_Dictionary, isa_proxy, 0,
                  sizeof(::RunConfig) );
      instance.SetNew(&new_RunConfig);
      instance.SetNewArray(&newArray_RunConfig);
      instance.SetDelete(&delete_RunConfig);
      instance.SetDeleteArray(&deleteArray_RunConfig);
      instance.SetDestructor(&destruct_RunConfig);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::RunConfig*)
   {
      return GenerateInitInstanceLocal((::RunConfig*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::RunConfig*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *RunConfig_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::RunConfig*)0x0)->GetClass();
      RunConfig_TClassManip(theClass);
   return theClass;
   }

   static void RunConfig_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_recocLcLMEMResult(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) ::reco::MEMResult : new ::reco::MEMResult;
   }
   static void *newArray_recocLcLMEMResult(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) ::reco::MEMResult[nElements] : new ::reco::MEMResult[nElements];
   }
   // Wrapper around operator delete
   static void delete_recocLcLMEMResult(void *p) {
      delete ((::reco::MEMResult*)p);
   }
   static void deleteArray_recocLcLMEMResult(void *p) {
      delete [] ((::reco::MEMResult*)p);
   }
   static void destruct_recocLcLMEMResult(void *p) {
      typedef ::reco::MEMResult current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::reco::MEMResult

namespace ROOT {
   // Wrappers around operator new
   static void *new_edmcLcLWrapperlEvectorlErecocLcLMEMResultgRsPgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) ::edm::Wrapper<vector<reco::MEMResult> > : new ::edm::Wrapper<vector<reco::MEMResult> >;
   }
   static void *newArray_edmcLcLWrapperlEvectorlErecocLcLMEMResultgRsPgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) ::edm::Wrapper<vector<reco::MEMResult> >[nElements] : new ::edm::Wrapper<vector<reco::MEMResult> >[nElements];
   }
   // Wrapper around operator delete
   static void delete_edmcLcLWrapperlEvectorlErecocLcLMEMResultgRsPgR(void *p) {
      delete ((::edm::Wrapper<vector<reco::MEMResult> >*)p);
   }
   static void deleteArray_edmcLcLWrapperlEvectorlErecocLcLMEMResultgRsPgR(void *p) {
      delete [] ((::edm::Wrapper<vector<reco::MEMResult> >*)p);
   }
   static void destruct_edmcLcLWrapperlEvectorlErecocLcLMEMResultgRsPgR(void *p) {
      typedef ::edm::Wrapper<vector<reco::MEMResult> > current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::edm::Wrapper<vector<reco::MEMResult> >

namespace ROOT {
   // Wrappers around operator new
   static void *new_edmcLcLReflEvectorlErecocLcLMEMResultgRcOrecocLcLMEMResultcOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMEMResultgRcOrecocLcLMEMResultgRsPgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) ::edm::Ref<vector<reco::MEMResult>,reco::MEMResult,edm::refhelper::FindUsingAdvance<vector<reco::MEMResult>,reco::MEMResult> > : new ::edm::Ref<vector<reco::MEMResult>,reco::MEMResult,edm::refhelper::FindUsingAdvance<vector<reco::MEMResult>,reco::MEMResult> >;
   }
   static void *newArray_edmcLcLReflEvectorlErecocLcLMEMResultgRcOrecocLcLMEMResultcOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMEMResultgRcOrecocLcLMEMResultgRsPgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) ::edm::Ref<vector<reco::MEMResult>,reco::MEMResult,edm::refhelper::FindUsingAdvance<vector<reco::MEMResult>,reco::MEMResult> >[nElements] : new ::edm::Ref<vector<reco::MEMResult>,reco::MEMResult,edm::refhelper::FindUsingAdvance<vector<reco::MEMResult>,reco::MEMResult> >[nElements];
   }
   // Wrapper around operator delete
   static void delete_edmcLcLReflEvectorlErecocLcLMEMResultgRcOrecocLcLMEMResultcOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMEMResultgRcOrecocLcLMEMResultgRsPgR(void *p) {
      delete ((::edm::Ref<vector<reco::MEMResult>,reco::MEMResult,edm::refhelper::FindUsingAdvance<vector<reco::MEMResult>,reco::MEMResult> >*)p);
   }
   static void deleteArray_edmcLcLReflEvectorlErecocLcLMEMResultgRcOrecocLcLMEMResultcOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMEMResultgRcOrecocLcLMEMResultgRsPgR(void *p) {
      delete [] ((::edm::Ref<vector<reco::MEMResult>,reco::MEMResult,edm::refhelper::FindUsingAdvance<vector<reco::MEMResult>,reco::MEMResult> >*)p);
   }
   static void destruct_edmcLcLReflEvectorlErecocLcLMEMResultgRcOrecocLcLMEMResultcOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMEMResultgRcOrecocLcLMEMResultgRsPgR(void *p) {
      typedef ::edm::Ref<vector<reco::MEMResult>,reco::MEMResult,edm::refhelper::FindUsingAdvance<vector<reco::MEMResult>,reco::MEMResult> > current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::edm::Ref<vector<reco::MEMResult>,reco::MEMResult,edm::refhelper::FindUsingAdvance<vector<reco::MEMResult>,reco::MEMResult> >

namespace ROOT {
   // Wrappers around operator new
   static void *new_edmcLcLRefProdlEvectorlErecocLcLMEMResultgRsPgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) ::edm::RefProd<vector<reco::MEMResult> > : new ::edm::RefProd<vector<reco::MEMResult> >;
   }
   static void *newArray_edmcLcLRefProdlEvectorlErecocLcLMEMResultgRsPgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) ::edm::RefProd<vector<reco::MEMResult> >[nElements] : new ::edm::RefProd<vector<reco::MEMResult> >[nElements];
   }
   // Wrapper around operator delete
   static void delete_edmcLcLRefProdlEvectorlErecocLcLMEMResultgRsPgR(void *p) {
      delete ((::edm::RefProd<vector<reco::MEMResult> >*)p);
   }
   static void deleteArray_edmcLcLRefProdlEvectorlErecocLcLMEMResultgRsPgR(void *p) {
      delete [] ((::edm::RefProd<vector<reco::MEMResult> >*)p);
   }
   static void destruct_edmcLcLRefProdlEvectorlErecocLcLMEMResultgRsPgR(void *p) {
      typedef ::edm::RefProd<vector<reco::MEMResult> > current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::edm::RefProd<vector<reco::MEMResult> >

namespace ROOT {
   // Wrappers around operator new
   static void *new_edmcLcLRefVectorlEvectorlErecocLcLMEMResultgRcOrecocLcLMEMResultcOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMEMResultgRcOrecocLcLMEMResultgRsPgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) ::edm::RefVector<vector<reco::MEMResult>,reco::MEMResult,edm::refhelper::FindUsingAdvance<vector<reco::MEMResult>,reco::MEMResult> > : new ::edm::RefVector<vector<reco::MEMResult>,reco::MEMResult,edm::refhelper::FindUsingAdvance<vector<reco::MEMResult>,reco::MEMResult> >;
   }
   static void *newArray_edmcLcLRefVectorlEvectorlErecocLcLMEMResultgRcOrecocLcLMEMResultcOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMEMResultgRcOrecocLcLMEMResultgRsPgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) ::edm::RefVector<vector<reco::MEMResult>,reco::MEMResult,edm::refhelper::FindUsingAdvance<vector<reco::MEMResult>,reco::MEMResult> >[nElements] : new ::edm::RefVector<vector<reco::MEMResult>,reco::MEMResult,edm::refhelper::FindUsingAdvance<vector<reco::MEMResult>,reco::MEMResult> >[nElements];
   }
   // Wrapper around operator delete
   static void delete_edmcLcLRefVectorlEvectorlErecocLcLMEMResultgRcOrecocLcLMEMResultcOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMEMResultgRcOrecocLcLMEMResultgRsPgR(void *p) {
      delete ((::edm::RefVector<vector<reco::MEMResult>,reco::MEMResult,edm::refhelper::FindUsingAdvance<vector<reco::MEMResult>,reco::MEMResult> >*)p);
   }
   static void deleteArray_edmcLcLRefVectorlEvectorlErecocLcLMEMResultgRcOrecocLcLMEMResultcOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMEMResultgRcOrecocLcLMEMResultgRsPgR(void *p) {
      delete [] ((::edm::RefVector<vector<reco::MEMResult>,reco::MEMResult,edm::refhelper::FindUsingAdvance<vector<reco::MEMResult>,reco::MEMResult> >*)p);
   }
   static void destruct_edmcLcLRefVectorlEvectorlErecocLcLMEMResultgRcOrecocLcLMEMResultcOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMEMResultgRcOrecocLcLMEMResultgRsPgR(void *p) {
      typedef ::edm::RefVector<vector<reco::MEMResult>,reco::MEMResult,edm::refhelper::FindUsingAdvance<vector<reco::MEMResult>,reco::MEMResult> > current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::edm::RefVector<vector<reco::MEMResult>,reco::MEMResult,edm::refhelper::FindUsingAdvance<vector<reco::MEMResult>,reco::MEMResult> >

namespace ROOT {
   // Wrappers around operator new
   static void *new_recocLcLIntegralResult(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) ::reco::IntegralResult : new ::reco::IntegralResult;
   }
   static void *newArray_recocLcLIntegralResult(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) ::reco::IntegralResult[nElements] : new ::reco::IntegralResult[nElements];
   }
   // Wrapper around operator delete
   static void delete_recocLcLIntegralResult(void *p) {
      delete ((::reco::IntegralResult*)p);
   }
   static void deleteArray_recocLcLIntegralResult(void *p) {
      delete [] ((::reco::IntegralResult*)p);
   }
   static void destruct_recocLcLIntegralResult(void *p) {
      typedef ::reco::IntegralResult current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::reco::IntegralResult

namespace ROOT {
   // Wrappers around operator new
   static void *new_RunConfig(void *p) {
      return  p ? new(p) ::RunConfig : new ::RunConfig;
   }
   static void *newArray_RunConfig(Long_t nElements, void *p) {
      return p ? new(p) ::RunConfig[nElements] : new ::RunConfig[nElements];
   }
   // Wrapper around operator delete
   static void delete_RunConfig(void *p) {
      delete ((::RunConfig*)p);
   }
   static void deleteArray_RunConfig(void *p) {
      delete [] ((::RunConfig*)p);
   }
   static void destruct_RunConfig(void *p) {
      typedef ::RunConfig current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::RunConfig

namespace ROOT {
   static TClass *vectorlErecocLcLMEMResultgR_Dictionary();
   static void vectorlErecocLcLMEMResultgR_TClassManip(TClass*);
   static void *new_vectorlErecocLcLMEMResultgR(void *p = 0);
   static void *newArray_vectorlErecocLcLMEMResultgR(Long_t size, void *p);
   static void delete_vectorlErecocLcLMEMResultgR(void *p);
   static void deleteArray_vectorlErecocLcLMEMResultgR(void *p);
   static void destruct_vectorlErecocLcLMEMResultgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<reco::MEMResult>*)
   {
      vector<reco::MEMResult> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<reco::MEMResult>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<reco::MEMResult>", -2, "vector", 214,
                  typeid(vector<reco::MEMResult>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlErecocLcLMEMResultgR_Dictionary, isa_proxy, 4,
                  sizeof(vector<reco::MEMResult>) );
      instance.SetNew(&new_vectorlErecocLcLMEMResultgR);
      instance.SetNewArray(&newArray_vectorlErecocLcLMEMResultgR);
      instance.SetDelete(&delete_vectorlErecocLcLMEMResultgR);
      instance.SetDeleteArray(&deleteArray_vectorlErecocLcLMEMResultgR);
      instance.SetDestructor(&destruct_vectorlErecocLcLMEMResultgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<reco::MEMResult> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const vector<reco::MEMResult>*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlErecocLcLMEMResultgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<reco::MEMResult>*)0x0)->GetClass();
      vectorlErecocLcLMEMResultgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlErecocLcLMEMResultgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlErecocLcLMEMResultgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<reco::MEMResult> : new vector<reco::MEMResult>;
   }
   static void *newArray_vectorlErecocLcLMEMResultgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<reco::MEMResult>[nElements] : new vector<reco::MEMResult>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlErecocLcLMEMResultgR(void *p) {
      delete ((vector<reco::MEMResult>*)p);
   }
   static void deleteArray_vectorlErecocLcLMEMResultgR(void *p) {
      delete [] ((vector<reco::MEMResult>*)p);
   }
   static void destruct_vectorlErecocLcLMEMResultgR(void *p) {
      typedef vector<reco::MEMResult> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<reco::MEMResult>

namespace {
  void TriggerDictionaryInitialization_MEMDataFormatsMEMDataFormats_xr_Impl() {
    static const char* headers[] = {
0    };
    static const char* includePaths[] = {
"/grid_mnt/vol__vol_U__u/llr/info/chiron_u/CMSSW_8_0_12/src",
"/cvmfs/cms.cern.ch/slc6_amd64_gcc530/cms/cmssw/CMSSW_8_0_12/src",
"/cvmfs/cms.cern.ch/slc6_amd64_gcc530/lcg/root/6.06.00-ikhhed4/include",
"/cvmfs/cms.cern.ch/slc6_amd64_gcc530/external/boost/1.57.0-ikhhed/include",
"/cvmfs/cms.cern.ch/slc6_amd64_gcc530/external/pcre/8.37/include",
"/cvmfs/cms.cern.ch/slc6_amd64_gcc530/external/bz2lib/1.0.6/include",
"/cvmfs/cms.cern.ch/slc6_amd64_gcc530/external/clhep/2.2.0.4-giojec/include",
"/cvmfs/cms.cern.ch/slc6_amd64_gcc530/external/gsl/1.16/include",
"/cvmfs/cms.cern.ch/slc6_amd64_gcc530/external/hepmc/2.06.07/include",
"/cvmfs/cms.cern.ch/slc6_amd64_gcc530/external/libuuid/2.22.2/include",
"/cvmfs/cms.cern.ch/slc6_amd64_gcc530/external/python/2.7.11-ikhhed/include/python2.7",
"/cvmfs/cms.cern.ch/slc6_amd64_gcc530/external/tbb/44_20151115oss/include",
"/cvmfs/cms.cern.ch/slc6_amd64_gcc530/external/xz/5.2.1/include",
"/cvmfs/cms.cern.ch/slc6_amd64_gcc530/external/zlib/1.2.8/include",
"/cvmfs/cms.cern.ch/slc6_amd64_gcc530/lcg/root/6.06.00-ikhhed4/include",
"/grid_mnt/vol__vol_U__u/llr/info/chiron_u/CMSSW_8_0_12/",
0
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "MEMDataFormatsMEMDataFormats_xr dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_Autoloading_Map;
namespace reco{class __attribute__((annotate("$clingAutoload$MEMDataFormats/MEMDataFormats/interface/MEMResult.h")))  MEMResult;}
namespace std{template <typename _Tp> class __attribute__((annotate("$clingAutoload$string")))  allocator;
}
namespace edm{namespace refhelper{template <typename C, typename T> struct __attribute__((annotate("$clingAutoload$MEMDataFormats/MEMDataFormats/interface/IntegralResult.h")))  FindUsingAdvance;
}}
namespace edm{namespace refhelper{template <typename REFV> struct __attribute__((annotate("$clingAutoload$MEMDataFormats/MEMDataFormats/interface/IntegralResult.h")))  FindRefVectorUsingAdvance;
}}
namespace edm{namespace refhelper{template <typename C, typename T> struct __attribute__((annotate("$clingAutoload$MEMDataFormats/MEMDataFormats/interface/IntegralResult.h")))  FindTrait;
}}
namespace edm{namespace refhelper{template <typename C> struct __attribute__((annotate("$clingAutoload$MEMDataFormats/MEMDataFormats/interface/IntegralResult.h")))  ValueTrait;
}}
namespace reco{class __attribute__((annotate("$clingAutoload$MEMDataFormats/MEMDataFormats/interface/IntegralResult.h")))  IntegralResult;}
class __attribute__((annotate("$clingAutoload$MEMDataFormats/MEMDataFormats/interface/RunConfig.h")))  RunConfig;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "MEMDataFormatsMEMDataFormats_xr dictionary payload"

#ifndef G__VECTOR_HAS_CLASS_ITERATOR
  #define G__VECTOR_HAS_CLASS_ITERATOR 1
#endif
#ifndef CMS_DICT_IMPL
  #define CMS_DICT_IMPL 1
#endif
#ifndef _REENTRANT
  #define _REENTRANT 1
#endif
#ifndef GNUSOURCE
  #define GNUSOURCE 1
#endif
#ifndef __STRICT_ANSI__
  #define __STRICT_ANSI__ 1
#endif
#ifndef GNU_GCC
  #define GNU_GCC 1
#endif
#ifndef _GNU_SOURCE
  #define _GNU_SOURCE 1
#endif
#ifndef CMSSW_GIT_HASH
  #define CMSSW_GIT_HASH "CMSSW_8_0_12"
#endif
#ifndef PROJECT_NAME
  #define PROJECT_NAME "CMSSW"
#endif
#ifndef PROJECT_VERSION
  #define PROJECT_VERSION "CMSSW_8_0_12"
#endif
#ifndef BOOST_SPIRIT_THREADSAFE
  #define BOOST_SPIRIT_THREADSAFE 1
#endif
#ifndef PHOENIX_THREADSAFE
  #define PHOENIX_THREADSAFE 1
#endif
#ifndef CMSSW_REFLEX_DICT
  #define CMSSW_REFLEX_DICT 1
#endif

#define _BACKWARD_BACKWARD_WARNING_H
#include "MEMDataFormats/MEMDataFormats/interface/IntegralResult.h"
#include "MEMDataFormats/MEMDataFormats/interface/MEMResult.h"
#include "MEMDataFormats/MEMDataFormats/interface/MEMResultFwd.h"
#include "MEMDataFormats/MEMDataFormats/interface/RunConfig.h"
#include "MEMDataFormats/MEMDataFormats/interface/PyRun2EventData_t.h"

#include "MEM/MEMProducer/plugins/MEM_nplet.h"
#include "DataFormats/Common/interface/Wrapper.h"

namespace { struct dictionary {
    reco::MEMResult      mr;
    reco::IntegralResult ir;
    RunConfig rc;
    
    reco::MEMResultCollection m1;
    edm::Wrapper<reco::MEMResultCollection> w1;
    edm::Ref<reco::MEMResultCollection> r1;
    edm::RefProd<reco::MEMResultCollection> rp1;
    edm::RefVector<reco::MEMResultCollection> rv1;
    
};
}
#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[]={
"RunConfig", payloadCode, "@",
"edm::Ref<std::vector<reco::MEMResult>,reco::MEMResult,edm::refhelper::FindUsingAdvance<std::vector<reco::MEMResult>,reco::MEMResult> >", payloadCode, "@",
"edm::Ref<vector<reco::MEMResult>,reco::MEMResult,edm::refhelper::FindUsingAdvance<vector<reco::MEMResult>,reco::MEMResult> >", payloadCode, "@",
"edm::RefProd<std::vector<reco::MEMResult> >", payloadCode, "@",
"edm::RefProd<vector<reco::MEMResult> >", payloadCode, "@",
"edm::RefVector<std::vector<reco::MEMResult>,reco::MEMResult,edm::refhelper::FindUsingAdvance<std::vector<reco::MEMResult>,reco::MEMResult> >", payloadCode, "@",
"edm::RefVector<vector<reco::MEMResult>,reco::MEMResult,edm::refhelper::FindUsingAdvance<vector<reco::MEMResult>,reco::MEMResult> >", payloadCode, "@",
"edm::Wrapper<std::vector<reco::MEMResult> >", payloadCode, "@",
"edm::Wrapper<vector<reco::MEMResult> >", payloadCode, "@",
"reco::IntegralResult", payloadCode, "@",
"reco::MEMResult", payloadCode, "@",
nullptr};

    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("MEMDataFormatsMEMDataFormats_xr",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_MEMDataFormatsMEMDataFormats_xr_Impl, {}, classesHeaders);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_MEMDataFormatsMEMDataFormats_xr_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_MEMDataFormatsMEMDataFormats_xr() {
  TriggerDictionaryInitialization_MEMDataFormatsMEMDataFormats_xr_Impl();
}
