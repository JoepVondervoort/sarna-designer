/* Browser-side geometry checks. Run through verify_workbench_browser.py --visual-audit. */
(() => {
  const scope=document.querySelector('dialog[open]') || document;
  const visible=node=>node.getClientRects().length && node.getBoundingClientRect().width>0 && node.getBoundingClientRect().height>0 && node.checkVisibility({checkVisibilityCSS:true});
  const identify=node=>node.id || node.getAttribute('aria-label') || node.textContent.trim().slice(0,70) || node.tagName;
  const violations=[], fields=[...scope.querySelectorAll('input:not([type=checkbox]):not([type=radio]):not([type=hidden]),textarea,select')].filter(visible);
  const original=document.activeElement;
  let labelChecks=0;
  for(const field of fields){
    const label=field.labels?.[0];if(!label || !visible(label))continue;
    const texts=[], walker=document.createTreeWalker(label,NodeFilter.SHOW_TEXT);
    while(walker.nextNode()){
      const text=walker.currentNode;
      if(!text.textContent.trim() || text.parentElement.closest('input,textarea,select'))continue;
      const range=document.createRange();range.selectNodeContents(text);
      texts.push(...[...range.getClientRects()].filter(r=>r.width>0 && r.height>0));
    }
    if(!texts.length)continue;
    field.focus({preventScroll:true});
    const r=field.getBoundingClientRect(), style=getComputedStyle(field);
    const ring=style.outlineStyle==='none'?0:Math.max(0,parseFloat(style.outlineWidth)+parseFloat(style.outlineOffset));
    const labelBottom=Math.max(...texts.map(t=>t.bottom));
    // Adjacent labels (e.g. compact controls) intentionally share the control's baseline.
    const horizontal=texts.some(t=>t.right>r.left+1 && t.left<r.right-1);
    if(horizontal){labelChecks++;if(r.top-ring<labelBottom+0.5)violations.push({kind:'label-focus-overlap',field:identify(field),gap:Math.round((r.top-ring-labelBottom)*10)/10});}
    for(let parent=field.parentElement;parent && parent!==document.body;parent=parent.parentElement){
      const css=getComputedStyle(parent);
      if(!['hidden','clip','auto','scroll'].includes(css.overflowX))continue;
      const box=parent.getBoundingClientRect(),left=box.left+parent.clientLeft,right=left+parent.clientWidth;
      if(r.left-ring<left-1 || r.right+ring>right+1)violations.push({kind:'clipped-focus-ring',field:identify(field),parent:identify(parent)});
      break;
    }
  }
  original?.focus?.({preventScroll:true});
  const maps=[...scope.querySelectorAll('#promoter-map,#gene-overview-map')].filter(visible);
  let svgLabels=0;
  for(const svg of maps){
    const labels=[...svg.querySelectorAll('text')].filter(visible).map(node=>({node,rect:node.getBoundingClientRect()}));
    svgLabels+=labels.length;
    const bounds=svg.getBoundingClientRect();
    for(let i=0;i<labels.length;i++){
      const a=labels[i];
      if(a.rect.left<bounds.left-1 || a.rect.right>bounds.right+1 || a.rect.top<bounds.top-1 || a.rect.bottom>bounds.bottom+1)violations.push({kind:'map-text-outside',map:svg.id,text:identify(a.node)});
      for(let j=i+1;j<labels.length;j++){
        const b=labels[j];
        const x=Math.min(a.rect.right,b.rect.right)-Math.max(a.rect.left,b.rect.left);
        const y=Math.min(a.rect.bottom,b.rect.bottom)-Math.max(a.rect.top,b.rect.top);
        if(x>1 && y>1)violations.push({kind:'map-text-overlap',map:svg.id,a:identify(a.node),b:identify(b.node)});
      }
    }
  }
  if(scope===document){
    const heading=document.querySelector('.workspace-panel:not([hidden]) > .workspace-heading h2,.workspace-panel:not([hidden]) #design-intro h1');
    const nav=document.getElementById('workflow-nav').getBoundingClientRect();
    if(heading && visible(heading)){
      const box=heading.getBoundingClientRect();
      if(box.bottom>nav.top && box.top<nav.bottom+4)violations.push({kind:'heading-under-sticky-navigation',heading:identify(heading)});
    }
  }
  if(document.documentElement.scrollWidth>innerWidth+1)violations.push({kind:'page-horizontal-overflow',width:innerWidth,content:document.documentElement.scrollWidth});
  return {fields:fields.length,labelChecks,svgLabels,violations};
})()
