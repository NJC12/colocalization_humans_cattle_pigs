// Execute the dashboard's real logic against a stub DOM and report which
// elements actually received content. Catches "function called but never
// defined", which node --check cannot see.
const fs=require("fs");
const html=fs.readFileSync(process.argv[2],"utf8");
const scripts=[...html.matchAll(/<script>([\s\S]*?)<\/script>/g)].map(m=>m[1]);
const data=scripts.filter(s=>s.trim().startsWith("window.PUB_"));
const logic=scripts.find(s=>s.includes("renderMatrix"));

const painted={};
function el(id){
  return { id,
    set innerHTML(v){ painted[id]=(v||"").length; }, get innerHTML(){ return ""; },
    set outerHTML(v){ painted[id]=(v||"").length; },
    querySelector:sel=>el(sel), querySelectorAll:()=>[],
    addEventListener(){}, closest:()=>null, dataset:{}, checked:false,
    getAttribute:()=>null, setAttribute(){}, textContent:"" };
}
global.window={ matchMedia:()=>({addEventListener(){}}) };
// querySelectorAll returns a real array: the page spreads it and reads .length,
// and the scroll-spy bails out cleanly on an empty one.
global.document={ getElementById:id=>el(id), querySelector:s=>el(s),
                  querySelectorAll:()=>[], addEventListener(){}, documentElement:{} };
global.IntersectionObserver=class{ constructor(){} observe(){} };
global.getComputedStyle=()=>({ getPropertyValue:()=> "#888888" });
global.matchMedia=global.window.matchMedia;

data.forEach(s=>eval(s));
try { eval(logic); }
catch(e){ console.log("RENDER THREW:", e.constructor.name+":", e.message); process.exit(1); }

const want=["prov","strip","#matrix thead","#matrix tbody","gridBox","gridStack","chartSpread","cbCats","cbArms","cbSubs","foot"];
let bad=0;
for(const k of want){
  const n=painted[k];
  const ok = n!==undefined && n>0;
  if(!ok) bad++;
  console.log(`  ${ok?"painted":"EMPTY  "}  ${k.padEnd(12)} ${n!==undefined?n+" chars":"never written"}`);
}
console.log(bad? `\n${bad} target(s) never received content` : "\nall targets painted");
process.exit(bad?1:0);
