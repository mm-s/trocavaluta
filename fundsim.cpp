#include <iostream>
#include <chrono>


using namespace std;

struct params {
	
	
};

struct user {
	double trade_rate{0}; //once a day
	double trade_amount{0};

	double daily_profit(const double& commtrade) const {
		return trade_amount*trade_rate+commtrade;
	}
};

struct webuser: user {
	webuser() {
	trade_rate=1; //once a day
	trade_amount=100;
	}
};
struct apiuser: user {
	apiuser() {
	trade_rate=20; //once a day
	trade_amount=50;
	}
};

//sell=1.36407
//buy=1.36428

struct spread {
	chrono::system_clock::time_point t;
	double bid;
	double ask;
};


#include <cex/engine.h>
#include <ctime>


struct forex: public vector<spread> {
	forex() {
		reserve(1000000);
		ifstream is("/home/marcos/fx/DAT_NT_EURGBP_T_LAST_201501_short.csv");
		while(!is.eof()) {
			string line;
			getline(is,line);
			istringstream ils(line);
			vector<string> tokens;
			string l;
			while(getline(ils,l,';')) {
				tokens.push_back(l);
			}
			if (tokens.size()!=3) continue;
			spread sp;
			std::tm tm;
			strptime(tokens[0].c_str(), "%Y%m%d %H%M%S", &tm);
			sp.t=std::chrono::system_clock::from_time_t(std::mktime(&tm));
			sp.bid=1/stod(tokens[1]);
			sp.ask=sp.bid-0.00050;
			emplace_back(sp);
		}
		//load
		i=begin();
	}
	void plot(std::ostream& os) const {
		for (auto& s:*this) {
			os << s.t.time_since_epoch().count() << " " << s.bid << endl;
		}

	}

	bool advance(int t) {
	}
	const spread& get_spread() const {
		return *i;
	}
	iterator i;
};
using namespace curex::cex;

struct manager {
	typedef curex::cex::string string;
	typedef curex::cex::ostream ostream;
	typedef curex::cex::id id;
	manager(engine& e, int cur1, double amount1, int cur2, double amount2): _e(e), _cur1(cur1), _cur2(cur2) {
		id nid{0};
		ostream nullos(0);
		_acpwd=L"8888888888";
		_acid=_e.new_account(_acpwd,nullos);
		_e.set_commission_discount(_acid,1); //100% discount for investor
		_e.deposit(_acid,_acpwd,cur1,amount1,nullos);
		_e.deposit(_acid,_acpwd,cur2,amount2,nullos);
		wcout << "investor acc id:" << _acid << endl;

	}
	manager(manager&&)=delete;
	manager(const manager&)=delete;
	~manager() {
		ostream nullos(0);
		_e.cancel_all_orders(_acid,_acpwd,nullos);
		wcout << "goodbye investor" << endl;
		_e.status(_acid,_acpwd,wcout);
	}


	struct omgr: curex::cex::engine::dorders {
		typedef curex::cex::efp efp;
		typedef map<id,double> orders;
		typedef map<efp,orders> queues;
		queues q;
		omgr(curex::cex::engine::dorders& b, int cid):curex::cex::engine::dorders(b) {
			for (auto& d:*this) {
				if (d.get_cid()!=cid) continue;
				if (d.key!=L'U') continue;
				auto i=q.find(d.get_encoded_rate());
				if (i==q.end()) {
					orders ords;
					ords.emplace(d._id,d.get_amount());
					q.emplace(d.get_encoded_rate(),move(ords));
					continue;
				}
				i->second.emplace(d._id,d.get_amount());
			}
		}
		void add(curex::cex::id id, auto rate, double amount) {
				//wcout << L"Adding " << id << L" " << rate << L" " << amount << endl;

				auto i=q.find(rate);
				if (i==q.end()) {
					orders ords;
					ords.emplace(id,amount);
					q.emplace(rate,move(ords));
					return;
				}
				i->second.emplace(id,amount);
		}
		void dump(curex::cex::ostream& os) const {
			for (auto ords:q) {
				wcout << ords.first << L" " << ords.second.size() << endl;
			}
		}
		bool exist_order(const efp& rate) const {
			auto i=q.find(rate);
			if (i==q.end()) return false;
			return !i->second.empty();
		}
		double volume(const efp& rate) const {
//	dump(wcout);
			auto i=q.find(rate);
			if (i==q.end()) return 0;
			double ans=0;
			for (auto& j:i->second) {
				ans+=j.second;
			}
			return ans;
		}
		double available() const {
			double ans=0;
			for (auto& i:q) {
			for (auto& j:i.second) {
				ans+=j.second;
			}
			}
			return ans;
		}
		void cancel_queue(engine& e, id acid, const string& pwd, auto rate) {
			ostream nullos(0);
			//wcout << L"Cancelling queue " << rate << endl;
			auto i=q.find(rate);
			if (i==q.end()) return;
			auto j=i->second.rbegin();
			while (j!=i->second.rend()) {
				e.cancel_order(acid,pwd,j->first,nullos);
				i->second.erase((++j).base());
				j=i->second.rbegin();
			}
			q.erase(i);
		}
		//returns amount not cancelled
		double cancel_amount(engine& e, id acid, const string& pwd, auto rate, double amount) {
			//wcout << L"Cancelling amount " << rate << L" " << amount << endl;
			ostream nullos(0);
			auto i=q.find(rate);
			if (i==q.end()) {
				dump(wcout);
				wcout << L"queue not found!! " << rate << endl;
				abort();				
				return amount;
			}
			auto j=i->second.rbegin();
			while (j!=i->second.rend()) {
				e.cancel_order(acid,pwd,j->first,nullos);
				amount-=j->second;
				i->second.erase((++j).base());
				if (amount<0.01) break;
				j=i->second.rbegin();
			}
			if (i->second.empty()) q.erase(i);
			return amount;
		}
	};

	struct plan: vector<pair<curex::cex::efp,double>> {
		typedef curex::cex::market::spots spots;
		typedef curex::cex::engine::tradeinput tradeinput;
		enum mode {
			bid, ask
		};
		plan(const auto& vr, const omgr& orders, int cid, double vol0_wallet, const auto& fx, mode mode_):_mode(mode_), /*_vr(vr), _spots(vr)*/ wallet(vol0_wallet) {
			double vol0_ords=orders.available();
			double vol0=vol0_wallet+vol0_ords;
			//wcout << L"$ avail wallet:" << vol0_wallet << L" orders:" << vol0_ords << L" total:" << vol0 << endl;
			//bid: front es el rate mas pequenio 
			//ask: front es el rate mas grande
			auto spot=mode_==ask?vr.lower_rate(fx.ask):vr.upper_rate(fx.bid);
//wcout << *vr.upper_rate(fx.ask) <<  L" " << fx.ask << L" " << fx.bid << L" " << *vr.lower_rate(fx.bid) << endl;
			double pvol=0.001;
			double vol=vol0;
			while(spot!=vr.end() && vol>0) {
				double am=vol0*pvol;
				if (am>vol) am=vol;
				vol-=am;
				emplace_back(pair<curex::cex::efp,double>(*spot,am));
				pvol*=1.5;
				if (pvol>1) pvol=1;
				mode_==ask?--spot:++spot;
			}
			//front values first
			if (mode_==ask) {
				sort(begin(),end(),[&](const auto& i,const auto& j) -> bool { return i>j; });
			}
			else {
				sort(begin(),end(),[&](const auto& i,const auto& j) -> bool{ return i<j; });
			}

		}
		
		mode _mode;	
		
		
		double wallet;
		//const spots& _spots;
		bool fund_wallet(engine&e,id acid, const string& pwd,omgr& om,double amount,int cid,auto limitrate) {
			//wcout << L"fund_wallet " << cid << L" " << amount << L" " << limitrate << endl;
			for (auto i=rbegin(); i!=rend(); ++i) {
				auto rate=i->first;
			//wcout << L"check rate " << rate << endl;
				if (_mode==ask) {
					if (rate>=limitrate) {
			//wcout << L"limit reached(bid) " << rate << " " << limitrate << endl;
						return false;
					}
				}
				else {
					if (rate<=limitrate) {
			//wcout << L"limit reached(ask) " << rate << " " << limitrate << endl;
						return false;
					}
				}
				amount=om.cancel_amount(e,acid,pwd,rate,amount); //returns amount left
				if (amount<0.1) {
					//wcout << L"left amount to fund " << amount << endl;
					return true;
				}
			}
			wcout << L"Should not be here " << size() << endl;
			return false;
		}
		bool exists(const curex::cex::efp& rate) const {
			for (auto i:*this)
				if (i.first==rate) return true;
			return false;
		}		
		void execute(engine&e, omgr& om, id acid, const string& pwd, int cid, const string& mkt) { //"AUD/GBP"
			vector<curex::cex::efp> queues_to_cancel;
			ostream nullos(0);

			for (const auto& pi:om.q) { //cancel all orders not in plan
				auto rate=pi.first;
				if (exists(rate)) continue;
				queues_to_cancel.push_back(rate);
			}			
			for (auto i:queues_to_cancel) {
				om.cancel_queue(e,acid,pwd,i);
			}

			for (const auto& pi:*this) { //from front to deep
				auto rate=pi.first;
				double expected_vol=pi.second;
			//om.dump(wcout);
				double cur_vol=om.volume(rate);
				//wcout << rate << L" expected " << expected_vol << L" cur " << cur_vol << endl;
				double diff=expected_vol-cur_vol;
				//wcout << L"diff " << diff << endl;
				while (diff<-1e-5) {
					om.cancel_amount(e,acid,pwd,rate,-diff);
					cur_vol=om.volume(rate);
					//wcout << rate << L" expected " << expected_vol << L" cur " << cur_vol << endl;
					diff=expected_vol-cur_vol;
					//wcout << L"diff " << diff << endl;
				}
				double wallet_funds=e.get_wallet(acid,pwd,cid,nullos);
				//wcout << rate << L" expected " << expected_vol << L" cur " << cur_vol << L" wallet " << wallet_funds << endl;
				if (diff>wallet_funds) {
					//wcout << L"fund wallet " << diff << L" " << wallet_funds << endl;
					if (!fund_wallet(e,acid,pwd,om,diff-wallet_funds,cid,rate)) abort();
					wallet_funds=e.get_wallet(acid,pwd,cid,nullos);
					//wcout << L" wallet " << wallet_funds << endl;
					assert(diff<wallet_funds);
				}
				if (diff>1e-5) {
					//e.liquidity(1,3, wcout);
					std::vector<tradeinput> input;
					{
					tradeinput ti;
					ti.code=L"S";
					ti.set_amount(diff);
					ti.set_cur(cid);
					ti.set_rate(rate);
					ti.rateunits=mkt;
					input.push_back(ti);	
					}
					curex::cex::ostringstream os;
					e.trade_low(acid, pwd, input, os);
					curex::cex::string line;
					curex::cex::istringstream is(os.str());
//wcout << os.str() << endl;
					//getline(is, line);					//U 1 S 500.00000 EUR 1.28755 EUR/GBP
					//wcout << "last line " << line << endl;
					curex::cex::string code;
					curex::cex::id id;
					curex::cex::string op;
					curex::cex::string am;
					is >> code;
					if(!code.empty()) {
						is >> id;
						is >> op;
						is >> am;
						if(code!=L"U") {
							wcout << os.str() << endl;
							abort();
						}
						assert(op==L"S");
	//wcout << am << L" " << curex::cex::decode(curex::cex::encode(diff)) << endl;
						assert(am==curex::cex::decode(curex::cex::encode(diff)));
						om.add(id,rate,diff);
					}
				}
			}
		}
		//const curex::cex::market::spots& _vr;
		void dump(const curex::cex::string& prefix, curex::cex::ostream& os) const {
			for (auto i:*this) {
				os << prefix << (_mode==ask?L"ask ":L"bid ") << curex::cex::decode(i.first) << L" " << i.second << endl;
			}
		}
	};

/*
*/

	void adjust_liquidity(const spread& fx) {
		//wcout << "spread: " << fx.bid << " " << fx.ask << endl;
		auto mkt=_e.get_market_str(eur::id,gbp::id);

		auto vr=_e.valid_rates(gbp::id,eur::id);
		vr.sort();
		//wcout << "VALID RATES" << endl;
		//vr.dump(wcout);
		ostream nullos(0);

		auto ws=_e.get_wallets(_acid,_acpwd,nullos);
		int mpsz=0;
		auto dords=_e.get_orders(_acid,_acpwd);
		{
		int cid=curex::cex::eur::id;
		//wcout << L"plan EUR" << endl;
		omgr ords(dords,cid);
	//ords.dump(wcout);
		double vol0_wallet=ws.find(cid)->second;
		plan p(vr,ords,cid,vol0_wallet,fx,plan::ask);
		//p.dump(L" ",wcout);
		if (mpsz<p.size()) mpsz=p.size();
		p.execute(_e,ords, _acid,_acpwd,cid,mkt); 

		}
		{
		int cid=curex::cex::gbp::id;
		//wcout << L"plan GBP" << endl;
		double vol0_wallet=ws.find(cid)->second;
		omgr ords(dords,cid);
	//ords.dump(wcout);
		plan p(vr,ords,cid,vol0_wallet,fx,plan::bid);
		//p.dump(L" ",wcout);
		if (mpsz<p.size()) mpsz=p.size();
		p.execute(_e,ords,_acid,_acpwd,cid,mkt); 
		}
		//_e.liquidity(eur::id,gbp::id,mpsz,wcout);
		//_e.status(_acid,_acpwd,wcout);
	}



	int _cur1, _cur2;
	engine& _e;
	id _acid;
	string _acpwd;
};
#include <mstd/algorithm>
struct trader {
	typedef curex::cex::string string;
	typedef curex::cex::id id;
	typedef curex::cex::efp efp;
	typedef curex::cex::ostream ostream;
	typedef curex::cex::engine::tradeinput tradeinput;
	typedef curex::cex::engine::updateinput updateinput;
	typedef curex::cex::ostringstream ostringstream;
	typedef curex::cex::istringstream istringstream;
	//cur1 weak, see acquire
	trader(engine& e, int cur1, double amount1, int cur2, double amount2, string mkt): _e(e), _cur1(cur1), _cur2(cur2), mkt(mkt), front_cur1(cur1,true), front_cur2(cur2,false), valid_rates(_e.valid_rates(cur1,cur2)) {
		_acpwd=L"8888888888";
		ostream nullos(0);

		_acid=_e.new_account(_acpwd,nullos);
		_e.deposit(_acid,_acpwd,cur1,amount1,nullos);
		_e.deposit(_acid,_acpwd,cur2,amount2,nullos);
		wcout << "trader acc id:" << _acid << endl;

		orders_cur1=new orders(_e,_acid,_acpwd,mkt,_cur1);
		orders_cur2=new orders(_e,_acid,_acpwd,mkt,_cur2);
		trades_cur1=new trades(_e,_acid,_acpwd,mkt,_cur1);
		trades_cur2=new trades(_e,_acid,_acpwd,mkt,_cur2);

		valid_rates.sort();

	}
	trader(trader&& other): _e(other._e), _cur1(other._cur1), _cur2(other._cur2), mkt(move(other.mkt)), front_cur1(move(other.front_cur1)), front_cur2(move(other.front_cur2)), orders_cur1(other.orders_cur1), orders_cur2(other.orders_cur2), trades_cur1(other.trades_cur1), trades_cur2(other.trades_cur2), _acpwd(move(other._acpwd)), _acid(other._acid), valid_rates(move(other.valid_rates)) {
		other.orders_cur1=0;
		other.orders_cur2=0;
		other.trades_cur1=0;
		other.trades_cur2=0;
	}
	trader(const trader& other)=delete;

	~trader() {
		if (orders_cur1!=0) {
			ostream nullos(0);
			_e.cancel_all_orders(_acid,_acpwd,nullos);
			wcout << "goodbye trader" << endl;
			_e.status(_acid,_acpwd,wcout);
		}
		//_e.status(_acid,_acpwd,wcout);
		delete orders_cur1;
		delete orders_cur2;
		delete trades_cur1;
		delete trades_cur2;
	}
	engine& _e;
	id _acid;
	string _acpwd;
	int _cur1, _cur2;
	string mkt;
	double ref_cur1,ref_cur2;

	curex::cex::market::spots valid_rates;

	double get_ref(int cid) const {
		if (cid==_cur1) return ref_cur1;
		return ref_cur2;
	}

	struct order {
		order(): _id(0), amount(0), rate(0) {}
		order(id id_,double am, efp rt): _id(id_), amount(am), rate(rt) {}
		id _id;
		double amount;
		efp rate;
		bool valid() const { return _id!=0; }
		void dump(ostream& os) const {
			os << _id << L" " << amount << L" " << rate << endl;
		} 
	};
	struct trade {
		trade(): amount(0), rate(0) {}
		trade(double am, efp rt): amount(am), rate(rt) {}
		double amount;
		efp rate;
		void dump(ostream& os) const {
			os << amount << L" " << rate << endl;
		} 
	};
	struct trades:std::vector<trade> {
		trades(engine& e, id acid, string pwd, string mkt, int cur): _e(e), acid(acid), pwd(pwd), _mkt(mkt), _cur(cur)  {

		}
		double get_sell_amount(const curex::cex::engine::dtrade& o) const {
			if (o.buy=='S') return o.get_amount();
			//B 100.0000 GBP 2.0000 AUD/GBP
			//buy cur in the numerator
			if (o.market.find(o.cur)==0) return o.get_amount()/o.get_rate();
			//buy cur in the denominator
			return o.get_amount()*o.get_rate();

		}
		void dump(ostream& os) const {
			for (auto& i: *this) i.dump(os);
		} 
		void update(const curex::cex::engine::dtrades& mo) {
			const auto& o=mo.filter_front(_cur);
			clear();
			for (const auto& i:o) {
				emplace_back(trade(get_sell_amount(i),i.get_encoded_rate()));
			}
		}
		//M 121 B 80.00000 AUD 1.38995 AUD/USD 57.55600 USD 8.00000 AUD
		engine& _e;
		id acid;
		string pwd;
		int _cur;
		string _mkt;
	};

	struct orders: std::vector<order> {

		orders(engine& e, id acid, string pwd, string mkt, int cur): _e(e), acid(acid), pwd(pwd), _mkt(mkt), _cur(cur)  {

		}
		double get_sell_amount(const curex::cex::engine::dorder& o) const {
			if (o.buy=='S') return o.get_amount();
			//B 100.0000 GBP 2.0000 AUD/GBP
			//buy cur in the numerator
			if (o.market.find(o.cur)==0) return o.get_amount()/o.get_rate();
			//buy cur in the denominator
			return o.get_amount()*o.get_rate();

		}
		void dump(ostream& os) const {
			for (auto& i: *this) i.dump(os);
		} 

		void update(const curex::cex::engine::dorders& mo) {
			const auto& o=mo.filter_front(_cur);
			clear();
			for (const auto& i:o) {
				emplace_back(order(i._id,get_sell_amount(i),i.get_encoded_rate()));
			}
		}
		double get_exposure() const {
			return accumulate(begin(),end(),0.0,[](const auto& a, const auto& b) { return a+b.amount; });		
		}
		void reduce_exposure(double amount) {
			if (size()!=1) {
				dump(wcout);
				assert(false);
			}
			const auto& o=*begin();
			updateinput ui;

			ui.set_order_id(o._id);
			ui.set_amount(o.amount-amount);
			ui.set_rate(o.rate);
			ostream nullos(0);
			_e.update(acid, pwd, vector<updateinput>{ui}, nullos);		
		}
		void sell(efp rate, double amount) {
			std::vector<tradeinput> input;
			{
			tradeinput ti;
			ti.code=L"S";
			ti.set_amount(amount);
			ti.set_cur(_cur);
			ti.set_rate(rate);
			ti.rateunits=_mkt;
			input.push_back(ti);	
			}
			//U 1002 S 100.00000 EUR 1.28740 EUR/GBP

			curex::cex::ostringstream os;
			_e.trade(acid, pwd, input, os);
			curex::cex::string line;
			curex::cex::istringstream is(os.str());
//wcout << os.str() << endl;
			//getline(is, line);					//U 1 S 500.00000 EUR 1.28755 EUR/GBP
			//wcout << "last line " << line << endl;
			curex::cex::string code;
			curex::cex::id id;
			curex::cex::string op;
			curex::cex::string am;
			is >> code;
			if(!code.empty()) {
				is >> id;
				is >> op;
				is >> am;
				assert(op==L"S");
				if(code==L"U") {
					push_back(order(id,curex::cex::decode_dbl(curex::cex::encode(am)),rate));
				}
				else if(code==L"M") {
					//update matches
				}
				
			}
		}
		void cancel_not_at(efp rate) {
			ostream nullos(0);

			mstd::filter(*this,[&](const auto&i) -> bool {
				if (i.rate==rate) return false;
				ostringstream os;
				_e.cancel_order(acid, pwd, i._id, os);
				istringstream is(os.str());
				curex::cex::string c;
				is >> c;
				if (c==L"C") return true;
				wcout << L"Unexpected result" << endl;
				abort();
				return false;
			});
		}
		engine& _e;
		id acid;
		string pwd;
		int _cur;
		string _mkt;
	};
	orders* orders_cur1;
	orders* orders_cur2;
	trades* trades_cur1;
	trades* trades_cur2;
	orders& get_orders(int cid) { 
		if (cid==_cur1) return *orders_cur1;
		return *orders_cur2;
	}
	trades& get_trades(int cid) { 
		if (cid==_cur1) return *trades_cur1;
		return *trades_cur2;
	}

	struct spot {
		spot(efp rate, double vol): rate(rate), vol(vol) {}
		efp rate;
		double vol;
	};
	struct front:vector<spot> {
		front(int cur, bool weak): cur(cur), weak( weak) {}
		front(front&&other): cur(other.cur), weak( other.weak), vector<spot>(other) {}
		bool weak;
		int cur;
		static spot nospot;
		const spot& dropback(double from) const {
			if (empty()) return nospot;
			return *begin();
		}
		void update(const curex::cex::engine::dliq& liq) {
			const auto& o=liq.filter_front(cur);
			clear();
			for (const auto& i:o) {
				emplace_back(spot(i.get_encoded_rate(),i.get_amount()));
			}
			//first entry is the closest to the gap
			sort(begin(),end(),[&](const auto&a,const auto&b) -> bool { return weak?a.rate>b.rate:a.rate<b.rate; });
		}
	};

	front front_cur1;
	front front_cur2;
	const front& get_front(int cid) const { 
		if (cid==_cur1) return front_cur1;
		return front_cur2;
	}

	efp get_sit_rate(int cid) const {
		 efp s=get_front(cid).dropback(get_ref(cid)).rate;
		 if (s!=0) return s;  //empty front
		 return get_front(cid).weak?*valid_rates.lower_rate(ref_cur1):*valid_rates.upper_rate(ref_cur2);
	}

	void risk(int cid, double target_exposure) {
		orders& o=get_orders(cid);
		efp sit=get_sit_rate(cid);
		
		o.cancel_not_at(sit);
		double currrent_exposure=o.get_exposure();
		double diff=target_exposure-currrent_exposure;
		if (diff<0.5) {
			o.reduce_exposure(-diff);
		}
		else if (diff>0.5) {
			o.sell(sit,diff);
		}
	}

	void risk() {
		double exposure_cur1=100;
		double exposure_cur2=100;
		risk(_cur1,exposure_cur1);
		risk(_cur2,exposure_cur2);
	}
	void hedge() {
		
	}
	void acquire(const spread& fx) {
		ref_cur1=fx.ask;
		ref_cur2=fx.bid;

		const curex::cex::engine::dorders dords=_e.get_orders(_acid,_acpwd);
		auto mdorders=dords.filter_market(mkt);
		orders_cur1->update(mdorders);
		orders_cur2->update(mdorders);

		const curex::cex::engine::dliq dliq=_e.get_liquidity(_cur1,_cur2,8);
		front_cur1.update(dliq);
		front_cur2.update(dliq);

	}
	void do_trade(const spread& fx) {
		acquire(fx);
		hedge();
		risk();
	}
};
trader::spot trader::front::nospot(0,0);

struct traders:vector<trader> {
	traders(engine& e, int num, int cur1, double amount1, int cur2, double amount2): _e(e) {
		auto mkt=e.get_market_str(cur1,cur2);
		reserve(num);
		for (int i=0; i<num; ++i)
			push_back(trader(_e,cur1,amount1,cur2,amount2, mkt));

	}
	engine& _e;

	void do_trade(const spread& fx) {

		for(auto& i:*this)
			i.do_trade(fx);
	}

};

void burnout(int argc, char** argv) {
	if (argc!=2) {
		cout << "fx" << endl;
		exit(1);
	}
	std::string cmd=argv[1];
	if (cmd=="fx") {
		forex fx;
		fx.plot(cout);
		exit(0);
	}

	typedef curex::cex::ostream ostream;
	//engine::home=;
	//engine::init();
	ostream nullos(0);

	engine e("",nullos);
	e.set_ruler(eur::id,gbp::id, 1.10,1.60,0.00005);
	forex fx;
	{
	manager mgr(e,eur::id,500000,gbp::id,500000);
	traders trs(e,1,eur::id,5000,gbp::id,5000);

	int ii=0;
	for (auto sp=fx.cbegin(); sp!=fx.cend(); ++sp) {
wcout << L"------------------------------------------------------" << ++ii << endl;
		//e.set_ref_rate(eur::id,gbp::id,sp->bid+(sp->ask-sp->bid)/2); //adjust edges of rates

		mgr.adjust_liquidity(*sp);

//		if (++ii==2) exit(0);
		trs.do_trade(*sp);
	}
	}
	
}


int main(int argc, char** argv) {

	burnout(argc,argv);
	return 0;
	double new_webusers_rate{10}; //10 a day
	double new_apiusers_rate{1}; //1 a day

	int active_webusers{0};
	int active_apiusers{0};
	
	double trade_commission=0.005;
	webuser wu;
	apiuser au;
	double profit_webusers{0};
	double profit_apiusers{0};

	for (int i=0; i<365; ++i) {
		active_webusers+=new_webusers_rate;
		active_apiusers+=new_apiusers_rate;
		profit_webusers+=wu.daily_profit(trade_commission);
		profit_apiusers+=au.daily_profit(trade_commission);
		double profit=profit_webusers+profit_apiusers;
		cout << "day " << i << " profit $" << profit << endl;
	}


return 0;
}
